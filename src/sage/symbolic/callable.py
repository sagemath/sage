"""
Callable Symbolic Expressions

EXAMPLES:

When you do arithmetic with::

    sage: f(x, y, z) = sin(x+y+z)
    sage: g(x, y) = y + 2*x
    sage: f + g
    (x, y, z) |--> 2*x + y + sin(x + y + z)

::

    sage: f(x, y, z) = sin(x+y+z)
    sage: g(w, t) = cos(w - t)
    sage: f + g
    (t, w, x, y, z) |--> cos(-t + w) + sin(x + y + z)

::

    sage: f(x, y, t) = y*(x^2-t)
    sage: g(x, y, w) = x + y - cos(w)
    sage: f*g
    (x, y, t, w) |--> (x^2 - t)*(x + y - cos(w))*y

::

    sage: f(x,y, t) = x+y
    sage: g(x, y, w) = w + t
    sage: f + g
    (x, y, t, w) |--> t + w + x + y

TESTS:

The arguments in the definition must be symbolic variables (:issue:`10747`)::

    sage: f(1)=2
    Traceback (most recent call last):
    ...
    SyntaxError: can...t assign to function call...

    sage: f(x,1)=2
    Traceback (most recent call last):
    ...
    SyntaxError: can...t assign to function call...

    sage: f(1,2)=3
    Traceback (most recent call last):
    ...
    SyntaxError: can...t assign to function call...

    sage: f(1,2)=x
    Traceback (most recent call last):
    ...
    SyntaxError: can...t assign to function call...

    sage: f(x,2)=x
    Traceback (most recent call last):
    ...
    SyntaxError: can...t assign to function call...
"""
import sage.rings.abc
from sage.symbolic.ring import SymbolicRing, SR
from sage.categories.pushout import ConstructionFunctor
from sage.structure.factory import UniqueFactory


######################################################################
#  Callable functions
######################################################################

class CallableSymbolicExpressionFunctor(ConstructionFunctor):
    def __init__(self, arguments):
        """
        A functor which produces a CallableSymbolicExpressionRing from
        the SymbolicRing.

        EXAMPLES::

            sage: from sage.symbolic.callable import CallableSymbolicExpressionFunctor
            sage: x,y = var('x,y')
            sage: f = CallableSymbolicExpressionFunctor((x,y)); f
            CallableSymbolicExpressionFunctor(x, y)
            sage: f(SR)
            Callable function ring with arguments (x, y)

            sage: loads(dumps(f))
            CallableSymbolicExpressionFunctor(x, y)
        """
        self._arguments = arguments
        from sage.categories.rings import Rings
        self.rank = 3
        ConstructionFunctor.__init__(self, Rings(), Rings())

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.callable import CallableSymbolicExpressionFunctor
            sage: x,y = var('x,y')
            sage: CallableSymbolicExpressionFunctor((x,y))
            CallableSymbolicExpressionFunctor(x, y)
        """
        return "CallableSymbolicExpressionFunctor%s" % repr(self.arguments())

    def merge(self, other):
        """
        EXAMPLES::

            sage: from sage.symbolic.callable import CallableSymbolicExpressionFunctor
            sage: x,y = var('x,y')
            sage: a = CallableSymbolicExpressionFunctor((x,))
            sage: b = CallableSymbolicExpressionFunctor((y,))
            sage: a.merge(b)
            CallableSymbolicExpressionFunctor(x, y)
        """
        arguments = self.unify_arguments(other)
        return CallableSymbolicExpressionFunctor(arguments)

    def __call__(self, R):
        """
        EXAMPLES::

            sage: from sage.symbolic.callable import CallableSymbolicExpressionFunctor
            sage: x,y = var('x,y')
            sage: a = CallableSymbolicExpressionFunctor((x,y))
            sage: a(SR)
            Callable function ring with arguments (x, y)
        """
        if R is not SR:
            raise ValueError("can only make callable symbolic expression rings from the Symbolic Ring")
        return CallableSymbolicExpressionRing(self.arguments())

    def arguments(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.callable import CallableSymbolicExpressionFunctor
            sage: x,y = var('x,y')
            sage: a = CallableSymbolicExpressionFunctor((x,y))
            sage: a.arguments()
            (x, y)
        """
        return self._arguments

    def unify_arguments(self, x):
        r"""
        Take the variable list from another ``CallableSymbolicExpression``
        object and compare it with the current ``CallableSymbolicExpression``
        object's variable list, combining them according to the following rules:

        Let ``a`` be ``self``'s variable list, let ``b`` be ``y``'s
        variable list.

        #. If ``a == b``, then the variable lists are
           identical, so return that variable list.

        #. If ``a`` `\neq` ``b``, then check if the first `n` items in
           ``a`` are the first `n` items in ``b``, or vice versa. If
           so, return a list with these `n` items, followed by the
           remaining items in ``a`` and ``b`` sorted together in
           alphabetical order.


        .. NOTE::

           When used for arithmetic between
           ``CallableSymbolicExpression``'s, these rules ensure that
           the set of ``CallableSymbolicExpression``'s will have
           certain properties. In particular, it ensures that the set
           is a *commutative* ring, i.e., the order of the input
           variables is the same no matter in which order arithmetic
           is done.

        INPUT:

        - ``x`` -- a ``CallableSymbolicExpression``

        OUTPUT: a tuple of variables

        EXAMPLES::

            sage: from sage.symbolic.callable import CallableSymbolicExpressionFunctor
            sage: x,y = var('x,y')
            sage: a = CallableSymbolicExpressionFunctor((x,))
            sage: b = CallableSymbolicExpressionFunctor((y,))
            sage: a.unify_arguments(b)
            (x, y)

        AUTHORS:

        - Bobby Moretti: thanks to William Stein for the rules
        """
        a = self.arguments()
        b = x.arguments()

        # Rule #1
        if [str(y) for y in a] == [str(z) for z in b]:
            return a

        # Rule #2
        new_list = []
        done = False
        i = 0
        while not done and i < min(len(a), len(b)):
            if repr(a[i]) == repr(b[i]):
                new_list.append(a[i])
                i += 1
            else:
                done = True

        temp = set()
        # Sorting remaining variables.
        for j in range(i, len(a)):
            if a[j] not in temp:
                temp.add(a[j])

        for j in range(i, len(b)):
            if b[j] not in temp:
                temp.add(b[j])

        new_list.extend(sorted(temp, key=repr))
        return tuple(new_list)


class CallableSymbolicExpressionRing_class(SymbolicRing, sage.rings.abc.CallableSymbolicExpressionRing):
    def __init__(self, arguments):
        """
        EXAMPLES:

        We verify that coercion works in the case where ``x`` is not an
        instance of SymbolicExpression, but its parent is still the
        SymbolicRing::

            sage: f(x) = 1
            sage: f*e
            x |--> e

        TESTS::

            sage: TestSuite(f.parent()).run(skip=['_test_divides'])
        """
        self._arguments = arguments
        SymbolicRing.__init__(self, SR)
        self._populate_coercion_lists_(coerce_list=[SR])
        self.symbols = SR.symbols  # Use the same list of symbols as SR

    def _coerce_map_from_(self, R):
        """
        EXAMPLES::

            sage: f(x,y) = x^2 + y
            sage: g(x,y,z) = x + y + z
            sage: f.parent().has_coerce_map_from(g.parent())
            False
            sage: g.parent().has_coerce_map_from(f.parent())
            True
        """
        if isinstance(R, CallableSymbolicExpressionRing_class):
            args = self.arguments()
            if all(a in args for a in R.arguments()):
                return True
            else:
                return False
        return SymbolicRing._coerce_map_from_(self, R)

    def construction(self):
        """
        EXAMPLES::

            sage: f(x,y) = x^2 + y
            sage: f.parent().construction()
            (CallableSymbolicExpressionFunctor(x, y), Symbolic Ring)
        """
        return (CallableSymbolicExpressionFunctor(self.arguments()), SR)

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: f(x) = x+1; g(y) = y+1
            sage: f.parent()(g)
            x |--> y + 1
            sage: g.parent()(f)
            y |--> x + 1
            sage: f(x) = x+2*y; g(y) = y+3*x
            sage: f.parent()(g)
            x |--> 3*x + y
            sage: g.parent()(f)
            y |--> x + 2*y
        """
        return SymbolicRing._element_constructor_(self, x)

    def _repr_(self):
        """
        String representation of ring of callable symbolic expressions.

        EXAMPLES::

            sage: R = CallableSymbolicExpressionRing(var('x,y,theta'))
            sage: R._repr_()
            'Callable function ring with arguments (x, y, theta)'

        We verify that :issue:`12298` has been fixed::

            sage: S = CallableSymbolicExpressionRing([var('z')])
            sage: S._repr_()
            'Callable function ring with argument z'
        """
        if len(self._arguments) == 0:
            return "Callable function ring with no named arguments"
        elif len(self._arguments) == 1:
            return "Callable function ring with argument {}".format(self._arguments[0])
        else:
            return "Callable function ring with arguments {}".format(self._arguments)

    def arguments(self):
        """
        Return the arguments of ``self``.

        The order that the variables appear in ``self.arguments()`` is
        the order that is used in evaluating the elements of ``self``.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: f(x,y) = 2*x+y
            sage: f.parent().arguments()
            (x, y)
            sage: f(y,x) = 2*x+y
            sage: f.parent().arguments()
            (y, x)
        """
        return self._arguments

    args = arguments

    def _repr_element_(self, x):
        """
        Return the string representation of the Expression ``x``.

        EXAMPLES::

            sage: f(y,x) = x + y
            sage: f
            (y, x) |--> x + y
            sage: f.parent()
            Callable function ring with arguments (y, x)
        """
        args = self.arguments()
        repr_x = SymbolicRing._repr_element_(self, x)
        if len(args) == 1:
            return "%s |--> %s" % (args[0], repr_x)
        else:
            args = ", ".join(map(str, args))
            return "(%s) |--> %s" % (args, repr_x)

    def _latex_element_(self, x):
        r"""
        Finds the LaTeX representation of this expression.

        EXAMPLES::

            sage: f(A, t, omega, psi) = A*cos(omega*t - psi)
            sage: f._latex_()
            '\\left( A, t, \\omega, \\psi \\right) \\ {\\mapsto} \\ A \\cos\\left(\\omega t - \\psi\\right)'

            sage: f(mu) =  mu^3
            sage: f._latex_()
            '\\mu \\ {\\mapsto}\\ \\mu^{3}'
        """
        from sage.misc.latex import latex
        args = self.args()
        args = [latex(arg) for arg in args]
        latex_x = SymbolicRing._latex_element_(self, x)
        if len(args) == 1:
            return r"%s \ {\mapsto}\ %s" % (args[0], latex_x)
        else:
            vars = ", ".join(args)
            return r"\left( %s \right) \ {\mapsto} \ %s" % (vars, latex_x)

    def _call_element_(self, _the_element, *args, **kwds):
        """
        Calling a callable symbolic expression returns a symbolic expression
        with the appropriate arguments substituted.

        EXAMPLES::

            sage: var('a, x, y, z')
            (a, x, y, z)
            sage: f(x,y) = a + 2*x + 3*y + z
            sage: f
            (x, y) |--> a + 2*x + 3*y + z
            sage: f(1,2)
            a + z + 8
            sage: f(y=2, a=-1)
            2*x + z + 5

        Note that keyword arguments will override the regular arguments.
        ::


            sage: f.arguments()
            (x, y)
            sage: f(1,2)
            a + z + 8
            sage: f(10,2)
            a + z + 26
            sage: f(10,2,x=1)
            a + z + 8
            sage: f(z=100)
            a + 2*x + 3*y + 100
        """
        if any(type(arg).__module__ == 'numpy' and type(arg).__name__ == "ndarray" for arg in args):  # avoid importing
            raise NotImplementedError("Numpy arrays are not supported as arguments for symbolic expressions")

        d = dict(zip([repr(_) for _ in self.arguments()], args))
        d.update(kwds)
        return SR(_the_element.substitute(**d))

    # __reduce__ gets replaced by the CallableSymbolicExpressionRingFactory
    __reduce__ = object.__reduce__


class CallableSymbolicExpressionRingFactory(UniqueFactory):
    def create_key(self, args, check=True):
        """
        EXAMPLES::

            sage: x,y = var('x,y')
            sage: CallableSymbolicExpressionRing.create_key((x,y))
            (x, y)
        """
        if check:
            from sage.structure.element import Expression
            if len(args) == 1 and isinstance(args[0], (list, tuple)):
                args, = args
            for arg in args:
                if not (isinstance(arg, Expression) and arg.is_symbol()):
                    raise TypeError("must construct a function with a tuple (or list) of variables")
            args = tuple(args)
        return args

    def create_object(self, version, key, **extra_args):
        """
        Return a CallableSymbolicExpressionRing given a version and a key.

        EXAMPLES::

            sage: x,y = var('x,y')
            sage: CallableSymbolicExpressionRing.create_object(0, (x, y))
            Callable function ring with arguments (x, y)
        """
        return CallableSymbolicExpressionRing_class(key)


CallableSymbolicExpressionRing = CallableSymbolicExpressionRingFactory('sage.symbolic.callable.CallableSymbolicExpressionRing')
