"""
Assumptions

The ``GenericDeclaration`` class provides assumptions about a symbol or
function in verbal form. Such assumptions can be made using the :func:`assume`
function in this module, which also can take any relation of symbolic
expressions as argument. Use :func:`forget` to clear all assumptions.
Creating a variable with a specific domain is equivalent with making an
assumption about it.

There is only rudimentary support for consistency and satisfiability checking
in Sage. Assumptions are used both in Maxima and Pynac to support or refine
some computations. In the following we show how to make and query assumptions.
Please see the respective modules for more practical examples.

In addition to the global :func:`assumptions` database, :func:`assuming`
creates reusable, stackable context managers allowing for temporary
updates of the database for evaluation of a (block of) statements.

EXAMPLES:

The default domain of a symbolic variable is the complex plane::

    sage: var('x')
    x
    sage: x.is_real()
    False
    sage: assume(x,'real')
    sage: x.is_real()
    True
    sage: forget()
    sage: x.is_real()
    False

Here is the list of acceptable features::

    sage: ", ".join(map(str, maxima("features")._sage_()))
    'integer, noninteger, even, odd, rational, irrational, real, imaginary,
    complex, analytic, increasing, decreasing, oddfun, evenfun, posfun,
    constant, commutative, lassociative, rassociative, symmetric,
    antisymmetric, integervalued'

Set positive domain using a relation::

    sage: assume(x>0)
    sage: x.is_positive()
    True
    sage: x.is_real()
    True
    sage: assumptions()
    [x > 0]

Assumptions also affect operations that do not use Maxima::

    sage: forget()
    sage: assume(x, 'even')
    sage: assume(x, 'real')
    sage: (-1)^x
    1
    sage: (-gamma(pi))^x
    gamma(pi)^x
    sage: binomial(2*x, x).is_integer()
    True

Assumptions are added and in some cases checked for consistency::

    sage: assume(x>0)
    sage: assume(x<0)
    Traceback (most recent call last):
    ...
    ValueError: Assumption is inconsistent
    sage: forget()
"""
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.rings.cc import CC
from sage.structure.element import Expression
from sage.structure.unique_representation import UniqueRepresentation

# #30074: We use the keys of a dict to store the assumptions.
# As of Python 3.6.x, dicts preserve the insertion order.
# In this way, we keep the same order of the assumptions
# as previous code that was using lists.
_assumptions = dict()

_valid_feature_strings = set()


class GenericDeclaration(UniqueRepresentation):
    """
    This class represents generic assumptions, such as a variable being
    an integer or a function being increasing. It passes such
    information to Maxima's declare (wrapped in a context so it is able
    to forget) and to Pynac.

    INPUT:

    - ``var`` -- the variable about which assumptions are being made

    - ``assumption`` -- string containing a Maxima feature, either user
      defined or in the list given by ``maxima('features')``

    EXAMPLES::

        sage: from sage.symbolic.assumptions import GenericDeclaration
        sage: decl = GenericDeclaration(x, 'integer')
        sage: decl.assume()
        sage: sin(x*pi)
        0
        sage: decl.forget()
        sage: sin(x*pi)
        sin(pi*x)
        sage: sin(x*pi).simplify()
        sin(pi*x)

    Here is the list of acceptable features::

        sage: ", ".join(map(str, maxima("features")._sage_()))
        'integer, noninteger, even, odd, rational, irrational, real, imaginary,
        complex, analytic, increasing, decreasing, oddfun, evenfun, posfun,
        constant, commutative, lassociative, rassociative, symmetric,
        antisymmetric, integervalued'

    Test unique representation behavior::

        sage: GenericDeclaration(x, 'integer') is GenericDeclaration(SR.var("x"), 'integer')
        True
    """

    def __init__(self, var, assumption):
        """
        This class represents generic assumptions, such as a variable being
        an integer or a function being increasing. It passes such
        information to maxima's declare (wrapped in a context so it is able
        to forget).

        INPUT:

        - ``var`` -- the variable about which assumptions are being made

        - ``assumption`` -- a Maxima feature, either user
          defined or in the list given by ``maxima('features')``

        EXAMPLES::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'integer')
            sage: decl.assume()
            sage: sin(x*pi)
            0
            sage: decl.forget()
            sage: sin(x*pi)
            sin(pi*x)

        Here is the list of acceptable features::

            sage: ", ".join(map(str, maxima("features")._sage_()))
            'integer, noninteger, even, odd, rational, irrational, real,
            imaginary, complex, analytic, increasing, decreasing, oddfun,
            evenfun, posfun, constant, commutative, lassociative, rassociative,
            symmetric, antisymmetric, integervalued'
        """
        self._var = var
        self._assumption = assumption
        self._context = None

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: GenericDeclaration(x, 'foo')
            x is foo
        """
        return "%s is %s" % (self._var, self._assumption)

    def has(self, arg):
        """
        Check if this assumption contains the argument ``arg``.

        EXAMPLES::

            sage: from sage.symbolic.assumptions import GenericDeclaration as GDecl
            sage: var('y')
            y
            sage: d = GDecl(x, 'integer')
            sage: d.has(x)
            True
            sage: d.has(y)
            False
        """
        return (arg - self._var).is_trivial_zero()

    def _validate_feature(self):
        """
        Check if this assumption is a known maxima feature, raise an error otherwise.

        EXAMPLES::

            sage: from sage.symbolic.assumptions import GenericDeclaration as GDecl
            sage: var('b')
            b
            sage: GDecl(b, 'bougie')
            b is bougie
            sage: _.assume()
            Traceback (most recent call last):
            ...
            ValueError: bougie not a valid assumption, must be one of ['analytic', ... 'symmetric']
        """
        from sage.calculus.calculus import maxima
        global _valid_feature_strings
        if self._assumption in _valid_feature_strings:
            return
        # We get the list here because features may be added with time.
        _valid_feature_strings.update(repr(x).strip() for x in list(maxima("features")))
        if self._assumption in _valid_feature_strings:
            return
        raise ValueError("%s not a valid assumption, must be one of %s"
                         % (self._assumption, sorted(_valid_feature_strings)))

    def assume(self):
        """
        Make this assumption.

        TESTS::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'even')
            sage: decl.assume()
            sage: cos(x*pi).simplify()
            1
            sage: decl2 = GenericDeclaration(x, 'odd')
            sage: decl2.assume()
            Traceback (most recent call last):
            ...
            ValueError: Assumption is inconsistent
            sage: decl.forget()
        """
        if self in _assumptions:
            return
        from sage.calculus.calculus import maxima
        cur = None
        context = None
        if self._context is None:
            self._validate_feature()
            cur = maxima.get("context")
            # newcontext makes a fresh context that only has $global as
            # a subcontext, and makes it the current $context,
            # but does not deactivate other current contexts.
            context = maxima.newcontext('context' + maxima._next_var_name())
            must_declare = True
        elif not maxima.featurep(self._var, self._assumption):
            # Reactivating a previously active context.
            # Run $declare again with the assumption
            # to catch possible inconsistency
            # with the active contexts.
            cur = maxima.get("context")
            # Redeclaring on the existing context does not seem to trigger
            # inconsistency checking.
            # maxima.set("context", self._context._maxima_init_())
            # Instead, use a temporary context for this purpose
            context = maxima.newcontext('context' + maxima._next_var_name())
            must_declare = True
        else:
            must_declare = False

        if must_declare:
            try:
                maxima.eval("declare(%s, %s)" % (self._var._maxima_init_(), self._assumption))
            except RuntimeError as mess:
                if 'inconsistent' in str(mess):  # note Maxima doesn't tell you if declarations are redundant
                    # Inconsistency with one of the active contexts.
                    raise ValueError("Assumption is inconsistent")
                else:
                    raise
            else:
                if self._context is None:
                    self._context = context
                    context = None
            finally:
                assert cur is not None
                maxima.set("context", cur)
                if context is not None:
                    maxima.killcontext(context)

        maxima.activate(self._context)
        self._var.decl_assume(self._assumption)
        _assumptions[self] = True

    def forget(self):
        """
        Forget this assumption.

        TESTS::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'odd')
            sage: decl.assume()
            sage: cos(pi*x)
            cos(pi*x)
            sage: cos(pi*x).simplify()
            -1
            sage: decl.forget()
            sage: cos(x*pi).simplify()
            cos(pi*x)
        """
        self._var.decl_forget(self._assumption)
        from sage.calculus.calculus import maxima
        if self._context is not None:
            try:
                del _assumptions[self]
            except KeyError:
                return
            maxima.deactivate(self._context)
        else:  # trying to forget a declaration explicitly rather than implicitly
            for x in _assumptions:
                if repr(self) == repr(x):  # so by implication x is also a GenericDeclaration
                    x.forget()
                    break
            return

    def contradicts(self, soln):
        """
        Return ``True`` if this assumption is violated by the given
        variable assignment(s).

        INPUT:

        - ``soln`` -- either a dictionary with variables as keys or a symbolic
          relation with a variable on the left hand side

        EXAMPLES::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: GenericDeclaration(x, 'integer').contradicts(x==4)
            False
            sage: GenericDeclaration(x, 'integer').contradicts(x==4.0)
            False
            sage: GenericDeclaration(x, 'integer').contradicts(x==4.5)
            True
            sage: GenericDeclaration(x, 'integer').contradicts(x==sqrt(17))
            True
            sage: GenericDeclaration(x, 'noninteger').contradicts(x==sqrt(17))
            False
            sage: GenericDeclaration(x, 'noninteger').contradicts(x==17)
            True
            sage: GenericDeclaration(x, 'even').contradicts(x==3)
            True
            sage: GenericDeclaration(x, 'complex').contradicts(x==3)
            False
            sage: GenericDeclaration(x, 'imaginary').contradicts(x==3)
            True
            sage: GenericDeclaration(x, 'imaginary').contradicts(x==I)
            False

            sage: var('y,z')
            (y, z)
            sage: GenericDeclaration(x, 'imaginary').contradicts(x==y+z)
            False

            sage: GenericDeclaration(x, 'rational').contradicts(y==pi)
            False
            sage: GenericDeclaration(x, 'rational').contradicts(x==pi)
            True
            sage: GenericDeclaration(x, 'irrational').contradicts(x!=pi)
            False
            sage: GenericDeclaration(x, 'rational').contradicts({x: pi, y: pi})
            True
            sage: GenericDeclaration(x, 'rational').contradicts({z: pi, y: pi})
            False
        """
        if isinstance(soln, dict):
            value = soln.get(self._var)
            if value is None:
                return False
        elif soln.lhs() == self._var:
            value = soln.rhs()
        else:
            return False
        try:
            CC(value)
        except TypeError:
            return False
        if self._assumption == 'integer':
            return value not in ZZ
        elif self._assumption == 'noninteger':
            return value in ZZ
        elif self._assumption == 'even':
            return value not in ZZ or bool(ZZ(value) % 2)
        elif self._assumption == 'odd':
            return value not in ZZ or not (ZZ(value) % 2)
        elif self._assumption == 'rational':
            return value not in QQ
        elif self._assumption == 'irrational':
            return value in QQ
        elif self._assumption == 'real':
            return value not in RR
        elif self._assumption == 'imaginary':
            return value not in CC or CC(value).real() != 0
        elif self._assumption == 'complex':
            return value not in CC


def preprocess_assumptions(args):
    """
    Turn a list of the form ``(var1, var2, ..., 'property')`` into a
    sequence of declarations ``(var1 is property), (var2 is property),
    ...``

    EXAMPLES::

        sage: from sage.symbolic.assumptions import preprocess_assumptions
        sage: preprocess_assumptions([x, 'integer', x > 4])
        [x is integer, x > 4]
        sage: var('x, y')
        (x, y)
        sage: preprocess_assumptions([x, y, 'integer', x > 4, y, 'even'])
        [x is integer, y is integer, x > 4, y is even]
    """
    args = list(args)
    last = None
    for i, x in reversed(list(enumerate(args))):
        if isinstance(x, str):
            del args[i]
            last = x
        elif ((not hasattr(x, 'assume')
               or (isinstance(x, Expression) and x.is_symbol()))
              and last is not None):
            args[i] = GenericDeclaration(x, last)
        else:
            last = None
    return args


def assume(*args):
    r"""
    Make the given assumptions.

    INPUT:

    - ``*args`` -- a variable-length sequence of assumptions, each
      consisting of:

      - any number of symbolic inequalities, like ``0 < x, x < 1``

      - a subsequence of variable names, followed by some property that
        should be assumed for those variables; for example, ``x, y, z,
        'integer'`` would assume that each of ``x``, ``y``, and ``z``
        are integer variables, and ``x, 'odd'`` would assume that ``x``
        is odd (as opposed to even).

      The two types can be combined, but a symbolic inequality cannot
      appear in the middle of a list of variables.

    OUTPUT: if everything goes as planned, there is no output

    If you assume something that is not one of the two forms above, then
    an :exc:`AttributeError` is raised as we try to call its ``assume``
    method.

    If you make inconsistent assumptions (for example, that ``x`` is
    both even and odd), then a :exc:`ValueError` is raised.

    .. WARNING::

        Do not use Python's chained comparison notation in assumptions.
        Python literally translates the expression ``0 < x < 1`` to
        ``(0 < x) and (x < 1)``, but the value of ``bool(0 < x)`` is
        ``False`` when ``x`` is a symbolic variable. Therefore, by the
        definition of Python's logical "and" operator, the entire expression
        is equal to ``0 < x``.

    EXAMPLES:

    Assumptions are typically used to ensure certain relations are
    evaluated as true that are not true in general.

    Here, we verify that for `x>0`, `\sqrt{x^2}=x`::

        sage: assume(x > 0)
        sage: bool(sqrt(x^2) == x)
        True

    This will be assumed in the current Sage session until forgotten::

        sage: bool(sqrt(x^2) == x)
        True
        sage: forget()
        sage: bool(sqrt(x^2) == x)
        False

    Another major use case is in taking certain integrals and limits
    where the answers may depend on some sign condition::

        sage: var('x, n')
        (x, n)
        sage: assume(n+1>0)
        sage: integral(x^n,x)
        x^(n + 1)/(n + 1)
        sage: forget()

    ::

        sage: var('q, a, k')
        (q, a, k)
        sage: assume(q > 1)
        sage: sum(a*q^k, k, 0, oo)
        Traceback (most recent call last):
        ...
        ValueError: Sum is divergent.
        sage: forget()
        sage: assume(abs(q) < 1)
        sage: sum(a*q^k, k, 0, oo)
        -a/(q - 1)
        sage: forget()

    An integer constraint::

        sage: n,P,r,r2 = SR.var('n, P, r, r2')
        sage: assume(n, 'integer')
        sage: c = P*e^(r*n)
        sage: d = P*(1+r2)^n
        sage: solve(c==d,r2)
        [r2 == e^r - 1]
        sage: forget()

    Simplifying certain well-known identities works as well::

        sage: n = SR.var('n')
        sage: assume(n, 'integer')
        sage: sin(n*pi)
        0
        sage: forget()
        sage: sin(n*pi).simplify()
        sin(pi*n)

    Instead of using chained comparison notation, each relationship
    should be passed as a separate assumption::

        sage: x = SR.var('x')
        sage: assume(0 < x, x < 1) # instead of assume(0 < x < 1)
        sage: assumptions()
        [0 < x, x < 1]
        sage: forget()

    If you make inconsistent or meaningless assumptions,
    Sage will let you know::

        sage: assume(x<0)
        sage: assume(x>0)
        Traceback (most recent call last):
        ...
        ValueError: Assumption is inconsistent
        sage: assume(x<1)
        Traceback (most recent call last):
        ...
        ValueError: Assumption is redundant
        sage: assumptions()
        [x < 0]
        sage: forget()
        sage: assume(x,'even')
        sage: assume(x,'odd')
        Traceback (most recent call last):
        ...
        ValueError: Assumption is inconsistent
        sage: forget()

    You can also use assumptions to evaluate simple
    truth values::

        sage: x, y, z = var('x, y, z')
        sage: assume(x>=y,y>=z,z>=x)
        sage: bool(x==z)
        True
        sage: bool(z<x)
        False
        sage: bool(z>y)
        False
        sage: bool(y==z)
        True
        sage: forget()
        sage: assume(x>=1,x<=1)
        sage: bool(x==1)
        True
        sage: bool(x>1)
        False
        sage: forget()

    TESTS:

    Test that you can do two non-relational
    declarations at once (fixing :issue:`7084`)::

        sage: var('m,n')
        (m, n)
        sage: assume(n, 'integer'); assume(m, 'integer')
        sage: sin(n*pi).simplify()
        0
        sage: sin(m*pi).simplify()
        0
        sage: forget()
        sage: sin(n*pi).simplify()
        sin(pi*n)
        sage: sin(m*pi).simplify()
        sin(pi*m)

    Check that positive integers can be created (:issue:`20132`)

        sage: x = SR.var('x', domain='positive')
        sage: assume(x, 'integer')
        sage: x.is_positive() and x.is_integer()
        True
        sage: forget()

        sage: x = SR.var('x', domain='integer')
        sage: assume(x > 0)
        sage: x.is_positive() and x.is_integer()
        True
        sage: forget()

        sage: assume(x, "integer")
        sage: assume(x > 0)
        sage: x.is_positive() and x.is_integer()
        True
        sage: forget()

    Ensure that an :exc:`AttributeError` is raised if we are given junk::

        sage: assume(3)
        Traceback (most recent call last):
        ...
        AttributeError: 'sage.rings.integer.Integer' object has no
        attribute 'assume'...

    Ensure that we can combine the two types of assumptions, as documented::

        sage: x,y = SR.var('x,y')
        sage: assume(x > 0, x, y, 'integer')
        sage: assumptions()
        [x > 0, x is integer, y is integer]
        sage: forget()
        sage: assume(x, y, 'integer', x > 0)
        sage: assumptions()
        [x is integer, y is integer, x > 0]
        sage: forget()

    Test that our WARNING block is accurate::

        sage: x = SR.var('x')
        sage: bool(0 < x)
        False
        sage: 0 < x < 1
        0 < x
        sage: assume(0 < x < 1)
        sage: assumptions()
        [0 < x]
        sage: forget()

    Check that :issue:`28538` is fixed::

        sage: x, y = SR.var('x, y')
        sage: assume(x > 0)
        sage: assume(y > 0)
        sage: bool(y*(x - y) == 0)
        False
    """
    for x in preprocess_assumptions(args):
        if isinstance(x, (tuple, list)):
            assume(*x)
        else:
            x.assume()


def forget(*args):
    """
    Forget the given assumption, or call with no arguments to forget
    all assumptions.

    Here an assumption is some sort of symbolic constraint.

    INPUT:

    - ``*args`` -- assumptions (default: forget all assumptions)

    EXAMPLES:

    We define and forget multiple assumptions::

        sage: forget()
        sage: var('x,y,z')
        (x, y, z)
        sage: assume(x>0, y>0, z == 1, y>0)
        sage: sorted(assumptions(), key=lambda x:str(x))
        [x > 0, y > 0, z == 1]
        sage: forget(x>0, z==1)
        sage: assumptions()
        [y > 0]
        sage: assume(y, 'even', z, 'complex')
        sage: assumptions()
        [y > 0, y is even, z is complex]
        sage: cos(y*pi).simplify()
        1
        sage: forget(y,'even')
        sage: cos(y*pi).simplify()
        cos(pi*y)
        sage: assumptions()
        [y > 0, z is complex]
        sage: forget()
        sage: assumptions()
        []
    """
    if len(args) == 0:
        _forget_all()
        return
    for x in preprocess_assumptions(args):
        if isinstance(x, (tuple, list)):
            forget(*x)
        else:
            try:
                x.forget()
            except KeyError:
                raise TypeError("forget not defined for objects of type '%s'" % type(x))


def assumptions(*args):
    """
    List all current symbolic assumptions.

    INPUT:

    - ``args`` -- list of variables which can be empty

    OUTPUT: list of assumptions on variables; if ``args`` is empty it returns all
      assumptions

    EXAMPLES::

        sage: var('x, y, z, w')
        (x, y, z, w)
        sage: forget()
        sage: assume(x^2+y^2 > 0)
        sage: assumptions()
        [x^2 + y^2 > 0]
        sage: forget(x^2+y^2 > 0)
        sage: assumptions()
        []
        sage: assume(x > y)
        sage: assume(z > w)
        sage: sorted(assumptions(), key=lambda x: str(x))
        [x > y, z > w]
        sage: forget()
        sage: assumptions()
        []

    It is also possible to query for assumptions on a variable independently::

        sage: x, y, z = var('x y z')
        sage: assume(x, 'integer')
        sage: assume(y > 0)
        sage: assume(y**2 + z**2 == 1)
        sage: assume(x < 0)
        sage: assumptions()
        [x is integer, y > 0, y^2 + z^2 == 1, x < 0]
        sage: assumptions(x)
        [x is integer, x < 0]
        sage: assumptions(x, y)
        [x is integer, x < 0, y > 0, y^2 + z^2 == 1]
        sage: assumptions(z)
        [y^2 + z^2 == 1]
    """
    if len(args) == 0:
        return list(_assumptions)

    result = []
    if len(args) == 1:
        result.extend([statement for statement in _assumptions
                       if statement.has(args[0])])
    else:
        for v in args:
            result += [statement for statement in list(_assumptions)
                       if str(v) in str(statement)]
    return result


def _forget_all():
    """
    Forget all symbolic assumptions.

    This is called by ``forget()``.

    EXAMPLES::

        sage: forget()
        sage: var('x,y')
        (x, y)
        sage: assume(x > 0, y < 0)
        sage: bool(x*y < 0)      # means definitely true
        True
        sage: bool(x*y > 0)      # might not be true
        False
        sage: forget()    # implicitly calls _forget_all
        sage: bool(x*y < 0)      # might not be true
        False
        sage: bool(x*y > 0)      # might not be true
        False

    TESTS:

    Check that :issue:`7315` is fixed::

        sage: var('m,n')
        (m, n)
        sage: assume(n, 'integer'); assume(m, 'integer')
        sage: sin(n*pi).simplify()
        0
        sage: sin(m*pi).simplify()
        0
        sage: forget()
        sage: sin(n*pi).simplify()
        sin(pi*n)
        sage: sin(m*pi).simplify()
        sin(pi*m)
    """
    global _assumptions
    if not _assumptions:
        return
    for x in list(_assumptions):
        # need to do this because x.forget() removes x from _assumptions
        x.forget()
    _assumptions = dict()


class assuming:
    """
    Temporarily modify assumptions.

    Create a context manager in which temporary assumptions are added
    (or substituted) to the current assumptions set.

    The set of possible assumptions and declarations  is the same as for
    :func:`assume`.

    This can be useful in interactive mode to discover the assumptions
    necessary to a given integration, or the exact solution to a system of
    equations.

    It can also be used to explore the branches of a :func:`cases()` expression.

    As with :func:`assume`, it is an error to add an assumption either redundant
    or inconsistent with the current assumption set (unless ``replace=True`` is
    used). See examples.

    INPUT:

    - ``*args`` -- assumptions (same format as for :func:`assume`)

    - ``replace`` -- boolean (default: ``False``); specifies whether the new
      assumptions are added to (default) or replace (if ``replace=True``) the
      current assumption set

    OUTPUT:

    A context manager useable in a ``with`` statement (see examples).

    EXAMPLES:

    Basic functionality : inside a :func:`with assuming:` block, Sage uses the
    updated assumptions database. After exit, the original database is
    restored. ::

        sage: var("x")
        x
        sage: forget(assumptions())
        sage: solve(x^2 == 4,x)
        [x == -2, x == 2]
        sage: with assuming(x > 0):
        ....:     solve(x^2 == 4,x)
        [x == 2]
        sage: assumptions()
        []

    The local assumptions can be stacked. We can use this functionality to
    discover incrementally the assumptions necessary to a given calculation
    (and by the way, to check that Sage's default integrator
    (Maxima's, that is), sometimes nitpicks for naught). ::

        sage: var("y,k,theta")
        (y, k, theta)
        sage: dgamma(y,k,theta)=y^(k-1)*e^(-y/theta)/(theta^k*gamma(k))
        sage: integrate(dgamma(y,k,theta),y,0,oo)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional constraints; using the 'assume' command before evaluation *may* help (example of legal syntax is 'assume(theta>0)', see `assume?` for more details)
        Is theta positive or negative?
        sage: a1=assuming(theta>0)
        sage: with a1:integrate(dgamma(y,k,theta),y,0,oo)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional constraints; using the 'assume' command before evaluation *may* help (example of legal syntax is 'assume(k>0)', see `assume?` for more details)
        Is k positive, negative or zero?
        sage: a2=assuming(k>0)
        sage: with a1,a2:integrate(dgamma(y,k,theta),y,0,oo)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional constraints; using the 'assume' command before evaluation *may* help (example of legal syntax is 'assume(k>0)', see `assume?` for more details)
        Is k an integer?
        sage: a3=assuming(k,"noninteger")
        sage: with a1,a2,a3:integrate(dgamma(y,k,theta),y,0,oo)
        1
        sage: a4=assuming(k,"integer")
        sage: with a1,a2,a4:integrate(dgamma(y,k,theta),y,0,oo)
        1

    As mentioned above, it is an error to try to introduce redundant or
    inconsistent assumptions. ::

        sage: assume(x > 0)
        sage: with assuming(x > -1): "I won't see this"
        Traceback (most recent call last):
        ...
        ValueError: Assumption is redundant

        sage: with assuming(x < -1): "I won't see this"
        Traceback (most recent call last):
        ...
        ValueError: Assumption is inconsistent
    """
    def __init__(self, *args, **kwds):
        r"""
        EXAMPLES::

            sage: forget()
            sage: foo = assuming(x>0)
            sage: foo.Ass
            (x > 0,)
            sage: bool(x>-1)
            False
        """
        self.replace = kwds.pop("replace", False)
        self.Ass = args

    def __enter__(self):
        r"""
        EXAMPLES::

            sage: forget()
            sage: foo = assuming(x>0)
            sage: bool(x>-1)
            False
            sage: foo.__enter__()
            sage: bool(x>-1)
            True
            sage: foo.__exit__()
            sage: bool(x>-1)
            False
        """
        if self.replace:
            self.OldAss = assumptions()
            forget(assumptions())
        assume(self.Ass)

    def __exit__(self, *args, **kwds):
        r"""
        EXAMPLES::

            sage: forget()
            sage: foo = assuming(x>0)
            sage: bool(x>-1)
            False
            sage: foo.__enter__()
            sage: bool(x>-1)
            True
            sage: foo.__exit__()
            sage: bool(x>-1)
            False
            sage: forget()
        """
        if self.replace:
            forget(assumptions())
            assume(self.OldAss)
        else:
            if len(self.Ass) > 0:
                forget(self.Ass)
