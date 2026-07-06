# sage.doctest: optional regina snappy
r"""
Abstract class for Python internal interfaces

This class contains common functionality of interfaces to packages that can be
installed (using pip) as a Python library (called a Python-CAS in the sequel)
such as Regina or SnapPy.

AUTHORS:

- Sebastian Oehms (2026): first version (refactored from regina.py)
"""

##############################################################################
#       Copyright (C) 2026 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.interfaces.interface import (
    Interface,
    InterfaceElement,
    InterfaceFunction,
    InterfaceFunctionElement,
)
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.misc.instancedoc import instancedoc


class PythonInternalInterface(ExtraTabCompletion, Interface):
    r"""
    Python internal interface.

    EXAMPLES::

        sage: K = Knots().from_table(8, 21)
        sage: Kr = regina(K); Kr
        <regina.Link: 8-crossing knot: ----++-- ( _5 _0 ^1 _2 _3 ^6 _7 ^3 _4 ^5 _6 ^7 ^0 _1 ^2 ^4 )>
        sage: Kr.knotSig()
        'iabcdbefcdghaefghRsgF+m'

    More examples can be found in the module header.
    """
    def __init__(self, name):
        r"""
        Python constructor.

        TESTS::

            sage: TestSuite(regina).run(skip=['_test_pickling', '_test_category'])
        """
        Interface.__init__(self, name)
        self._initialized = False  # done lazily
        # Namespace of the partner Python project
        self._namespace = None
        # modules that contain class declarations of the partner Python project
        self._interface_modules = []
        # global variables of the interface including the namespace
        self._interface_globals = {}

    def _lazy_init(self):
        r"""
        Initialize the Python-CAS interpreter.

        Implemented according to R interface.

        EXAMPLES::

            sage: snappy._lazy_init()
        """
        if not self._initialized:
            self._initialized = True
            self._start()

    def __reduce__(self):
        r"""
        Helper for pickling.

        EXAMPLES::

            sage: p = regina.Polynomial([-2, ~7])  # indirect doctest
            sage: loads(dumps(p)) == p
            True
        """
        return self.__class__, ()

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: regina.an_element()      # indirect doctest
            2
        """
        return self(2)

    def _start(self):
        """
        Start up the interpreter and sets the initial prompt and options.

        This is called the first time the interface is actually used.

        EXAMPLES::

            sage: regina._start()
            sage: regina._namespace.Cyclotomic
            <class 'regina.engine.Cyclotomic'>
            sage: snappy._start()
            sage: snappy._namespace.Manifold
            <class 'SnapPy.Manifold'>
        """
        raise NotImplementedError('this method must be implemented in the child class')

    def _install_hints(self):
        """
        Hints for installing the interface on your computer.

        EXAMPLES::

            sage: len(regina._install_hints())
            99
        """
        raise NotImplementedError('this method must be implemented in the child class')

    def _eval(self, code):
        """
        Evaluates a command inside the Python-CAS interpreter and returns the output
        as a Python-CAS result.

        EXAMPLES::

            sage: regina._eval("Link('iabcdbefcdghaefghRsgF+m')")
            <regina.Link: 8-crossing knot: ++++--++ ( ^0 ^1 _2 ^3 _1 ^4 _5 ^2 _3 _6 ^7 _0 _4 ^5 ^6 _7 )>
            sage: snappy._eval("Link('9_42')")
            <Link 9_42: 1 comp; 9 cross>
        """
        self._lazy_init()
        globs = self._interface_globals
        if type(code) is str and code.find('=') < 0:
            pre_bracket = code.split('(')[0]
            if pre_bracket in globs:
                val = globs[pre_bracket]
                nam = self._namespace.__name__
                if callable(val):
                    return eval('%s.%s' % (nam, code), globs)
                return val
            return eval(code, globs)
        return exec(code, globs)

    def eval(self, code, *args, **kwds):
        """
        Evaluates a command inside the Python-CAS interpreter and returns the output
        in printable form.

        EXAMPLES::

            sage: regina.eval('1+1')
            '2'
        """
        return repr(self._eval(code))

    def get(self, var):
        """
        Get the value of the variable ``var``.

        EXAMPLES::

            sage: regina.get('Link')
            <class 'regina.engine.Link'>
            sage: snappy.get('Triangulation')
            <class 'SnapPy.Triangulation'>
        """
        return self._interface_globals[var]

    def set(self, var, value):
        """
        Set the variable ``var`` to the given ``value``.

        EXAMPLES::

            sage: regina.set('myLink', 'Link')
            sage: regina.get('myLink')
            <class 'regina.engine.Link'>
            sage: snappy.set('K9_15', 'Link("9_15")')
            sage: snappy.get('K9_15')
            <Link 9_15: 1 comp; 9 cross>
        """
        self._lazy_init()
        globs = self._interface_globals
        if not isinstance(value, str):
            globs[var] = value
            return
        try:
            val = self._eval(value)
            globs[var] = val
        except (NameError, AttributeError, KeyError):
            pass
        super().set(var, value)

    def _internal_namespace_object(self, x) -> bool:
        r"""
        Return ``True`` if ``x`` is an instance of a class
        from Python-CAS's namespace.

        EXAMPLES::

            sage: L = regina.get('Link')
            sage: regina._internal_namespace_object(L)
            False
            sage: regina._internal_namespace_object(L())
            True
            sage: snappy._internal_namespace_object([])
            False
        """
        self._lazy_init()
        ns = self._namespace.__dict__
        cl = x.__class__
        cln = cl.__name__
        if cln in ns:
            if cl == ns[cln]:
                return True
        im = self._interface_modules
        clm = cl.__module__.split('.')[0]
        if clm in im:
            return True
        return False

    def _coerce_from_special_method(self, x):
        """
        Try to coerce to ``self`` by calling a special underscore method.

        This method is overloaded to record the Sage parent of ``x`` in the
        interface element.

        EXAMPLES::

            sage: R.<u, v> = LaurentPolynomialRing(ZZ)
            sage: p = u*~v^3 + 3*v*~u + 5*u - 7
            sage: rp = regina(p)   # indirect doctest
            sage: rp._sage_parent
            Multivariate Laurent Polynomial Ring in u, v over Integer Ring
        """
        res = super()._coerce_from_special_method(x)
        if hasattr(x, 'parent'):
            res._sage_parent = x.parent()
        return res

    def _coerce_impl(self, x, use_special=True):
        r"""
        Coerce pure Python types via corresponding Sage objects.

        This method is overloaded to add Python types from self._namespace

        EXAMPLES::

            sage: L = regina.get('Link')
            sage: regina._coerce_impl(L())
            <regina.Link: Empty link>
        """
        if self._internal_namespace_object(x):
            return self(self._create(x))
        return super()._coerce_impl(x, use_special=use_special)

    def _convert_args_kwds(self, *args, **kwds):
        """
        Convert all of the ``args`` and ``kwds`` to instances of Python-CAS
        classes.

        EXAMPLES::

            sage: a = [regina(i) for i in range(3)]
            sage: b = list(range(3))
            sage: L = regina.get('Link')
            sage: C = L.fromKnotSig('iabcdbefcdghaefghRsgF+m')
            sage: D = (3,7)
            sage: regina._convert_args_kwds(a, b, C=C, D=D)
            (([0, 1, 2], [0, 1, 2]),
            {'C': <regina.Link: 8-crossing knot: ++++--++ ( ^0 ^1 _2 ^3 _1 ^4 _5 ^2 _3 _6 ^7 _0 _4 ^5 ^6 _7 )>,
             'D': (3, 7)})
        """
        def convert_arg(arg):
            coerce_name = '_%s_' % self.name()
            if isinstance(arg, InterfaceElement) and arg.parent() is self:
                return arg._inst
            if isinstance(arg, (list, tuple)):
                return type(arg)([convert_arg(i) for i in arg])
            if hasattr(arg, coerce_name):
                coerce = arg.__getattribute__(coerce_name)
                reg = coerce(self)
                return convert_arg(reg)
            return arg

        if args:
            args = list(args)
            for i, arg in enumerate(args):
                args[i] = convert_arg(arg)
        if kwds:
            for key, value in kwds.items():
                kwds[key] = convert_arg(value)
        return tuple(args), kwds

    def _function_call(self, name, *args, **kwds):
        r"""
        Perform a function call.

        EXAMPLES::

            sage: regina._function_call(regina.Polynomial, (-3, 5/3))
            <regina.PolynomialRational: 5/3 x - 3>
        """
        args, kwds = self._convert_args_kwds(*args, **kwds)
        if len(args) == 0:
            if len(kwds) == 0:
                res = name()
            else:
                res = name(**kwds)
        elif len(args) == 1:
            if len(kwds) == 0:
                res = name(args[0])
            else:
                res = name(args[0], **kwds)
        elif len(kwds) == 0:
            res = name(*args)
        else:
            res = name(*args, **kwds)

        # read back new values of the arguments and keywords
        def read_back(arg):
            if isinstance(arg, self._object_class()):
                self.set(arg._name, arg._inst)
        for arg in args:
            read_back(arg)
        for val in kwds.values():
            read_back(val)
        if res is not None:
            if self._internal_namespace_object(res):
                from sage.interfaces.regina import Regina
                if isinstance(self, Regina):
                    return self(res.__class__(res))  # this is the way to get a copy of a Regina object
                return self._object_class()(self, res)
            if type(res) in (list, tuple):
                return self._object_class()(self, res)
            return res

    def _equality_symbol(self):
        r"""
        EXAMPLES::

            sage: regina._equality_symbol()
            '=='
        """
        return '=='

    def _object_class(self):
        r"""
        Return the element class of this parent.
        This is used in the interface class.

        EXAMPLES::

            sage: regina._object_class()
            <class 'sage.interfaces.regina.ReginaElement'>
            sage: snappy._object_class()
            <class 'sage.interfaces.snappy.SnapPyElement'>
        """
        return PythonInternalElement

    def help(self, cmd, long=False):
        r"""
        Return the documentation of the given command via the Python internal
        interface.

        EXAMPLES::

            sage: regina.help('AbelianGroup')
            Represents a finitely generated abelian group.
            <BLANKLINE>
            The torsion elements of the group are stored in terms of their
            invariant factors. For instance, Z_2+Z_3 will appear as Z_6, and
            Z_2+Z_2+Z_3 will appear as Z_2+Z_6.
            <BLANKLINE>
            In general the factors will appear as Z_*d0*+...+Z_*dn*, where the
            invariant factors *di* are all greater than 1 and satisfy
            *d0*|*d1*|...|*dn*. Note that this representation is unique.
            <BLANKLINE>
            This class implements C++ move semantics and adheres to the C++
            Swappable requirement. It is designed to avoid deep copies wherever
            possible, even when passing or returning objects by value.

            sage: snappy.help('AbelianGroup')
            <BLANKLINE>
            An AbelianGroup object represents a finitely generated abelian group,
            usually the first homology group of a snappy Manifold.
            <BLANKLINE>
            Instantiate an abelian group by its elementary divisors:
            ...
        """
        self._lazy_init()
        dic = self._namespace.__dict__
        if cmd in dic:
            return print(dic[cmd].__doc__)
        raise NotImplementedError('no documentation available for %s' % cmd)

    def _tab_completion(self):
        r"""
        Return a list of all classes available through the interface.

        .. NOTE::

           Currently returns all keys of the namespace dictionary.

        EXAMPLES::

            sage: 'AbelianGroup' in regina._tab_completion()
            True
        """
        self._lazy_init()
        return list(self._namespace.__dict__)

    def __getattr__(self, attrname):
        r"""

        EXAMPLES::

            sage: type(regina.AbelianGroup)
            <class 'sage.interfaces.python_internal.PythonInternalFunction'>
            sage: regina.AbelianGroup._name
            <class 'regina.engine.AbelianGroup'>
        """
        if attrname[:1] == "_":
            raise AttributeError
        self._lazy_init()
        try:
            attr = self._namespace.__dict__[attrname]
        except KeyError:
            raise AttributeError
        if callable(attr):
            return PythonInternalFunction(self, self._namespace.__dict__[attrname])
        res = self(self._create(attr))
        res._inst = attr
        return res


@instancedoc
class PythonInternalElement(ExtraTabCompletion, InterfaceElement):
    r"""
    Element class of the Python internal interface.

    Its instances are usually constructed via the instance call of its parent.
    It wrapes the Python internal library for this object. In a session Python
    internal methods can be obtained using tab completion.

    EXAMPLES::

        sage: b = BraidGroup(3)((1,2,-1))
        sage: re = regina(b); re
        <regina.GroupExpression: g0 g1 g0^-1>
        sage: type(re)
        <class 'sage.interfaces.regina.ReginaElement'>
        sage: P = re.parent(); P
        Regina
        sage: type(P)
        <class 'sage.interfaces.regina.Regina'>

    Access to the Python-CAS expression objects::

        sage: res = re._inst
        sage: type(res)
         <class 'regina.engine.GroupExpression'>

    Applying Python-CAS methods::

        sage: re.cycleLeft(); re
        <regina.GroupExpression: g0^-1 g0 g1>

    Conversion to Sage::

        sage: re.sage() == b
        False
        sage: re.cycleRight()
        sage: re.sage() == b
        True

    TESTS::

        sage: p = regina.Polynomial([-2, ~7])
        sage: TestSuite(p).run(skip='_test_category')
    """

    _sage_parent = None  # for interface elements create from sage objects

    def _tab_completion(self):
        r"""
        Return a list of all methods of this object.

        EXAMPLES::

            sage: 'detail' in regina.AbelianGroup()._tab_completion()
            True
        """
        return dir(self._inst)

    def __getitem__(self, n):
        r"""
        EXAMPLES::

            sage: t = regina([-2, 7/3])
            sage: type(t)
            <class 'sage.interfaces.regina.ReginaElement'>
            sage: t[1]     # indirect doctest
            7/3
            sage: type(t[1])
            <class 'sage.interfaces.regina.ReginaElement'>
            sage: type(t[1]._inst)
            <class 'regina.engine.Rational'>
        """
        P = self.parent()
        return P._object_class()(P, self._inst[n])

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: list(regina([-2, 7/3]))
            [-2, 7/3]
        """
        for i in range(len(self)):
            yield self[i]

    def __getattr__(self, attrname):
        r"""
        EXAMPLES::

            sage: type(regina.AbelianGroup().detail)
            <class 'sage.interfaces.python_internal.PythonInternalFunctionElement'>
            sage: snappy.Link().PD_code
            <bound method Link.PD_code of <Link: 0 comp; 0 cross>>
        """
        P = self._check_valid()
        if attrname == '_inst':
            self._inst = P.get(self.name())
            return self._inst
        if attrname[:1] == "_":
            raise AttributeError
        else:
            inst = self._inst
            if hasattr(inst, attrname):
                attr = inst.__getattribute__(attrname)
                if callable(attr):
                    return PythonInternalFunctionElement(self, attr)
        raise AttributeError

    def __bool__(self):
        """
        Return whether this element is not ``False``.

        EXAMPLES::

            sage: M = regina.Matrix2(); M
            <regina.Matrix2: [[ 0 0 ] [ 0 0 ]]>
            sage: M.isIdentity()
            False
            sage: M.isZero()
            True
        """
        return bool(self._inst)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: a = regina(~5)
            sage: latex(a)
            \frac{1}{5}
            sage: b = regina(BraidGroup(4)((1, 2, 3, -2, -1))); b
            <regina.GroupExpression: g0 g1 g2 g1^-1 g0^-1>
            sage: latex(b)
            g_{0}g_{1}g_{2}g_{1}^{-1}g_{0}^{-1}
        """
        if hasattr(self, 'tex'):
            return self.tex()
        return super()._latex_()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: Q = regina('GroupPresentation()')
            sage: repr(Q)
            '<regina.GroupPresentation: < >>'
        """
        return repr(self._inst)

    def _richcmp_(self, other, op):
        r"""
        Comparison of interface elements.

        EXAMPLES::

            sage: a = regina(3)
            sage: b = regina(5)
            sage: a < b
            True
            sage: a > b
            False
            sage: F3 = FreeGroup(3)
            sage: a, b, c = F3.gens()
            sage: f = a**2*b*~c
            sage: fr = regina(f)
            sage: cr = regina(c)
            sage: fr == cr
            False
            sage: fr != fr * cr
            True
        """
        from sage.structure.richcmp import op_EQ, op_NE, rich_to_bool
        if self._inst == other._inst:
            return rich_to_bool(op, 0)
        if op == op_EQ:
            return False
        if op == op_NE:
            return True
        try:
            if self._inst < other._inst:
                return rich_to_bool(op, -1)
            if self._inst > other._inst:
                return rich_to_bool(op, 1)
        except TypeError:
            pass
        return super()._richcmp_(other, op)

    def _operation(self, operation, other=None):
        r"""
        Return the result of applying the binary operation
        ``operation`` on the arguments ``self`` and ``other``, or the
        unary operation on ``self`` if ``other`` is not given.

        This is a utility function which factors out much of the
        commonality used in the arithmetic operations for interface
        elements.

        INPUT:

        - ``operation`` -- string representing the operation
          being performed; for example, '*', or '1/'

        - ``other`` -- the other operand; if ``other`` is ``None``,
          then the operation is assumed to be unary rather than binary

        OUTPUT: an interface element

        EXAMPLES::

            sage: l = regina(range(2, 5))
            sage: [type(i) for i in l]
            [<class 'sage.interfaces.regina.ReginaElement'>,
             <class 'sage.interfaces.regina.ReginaElement'>,
             <class 'sage.interfaces.regina.ReginaElement'>]
            sage: s = sum(l); s, type(s)
            (9, <class 'sage.interfaces.regina.ReginaElement'>)
            sage: p = prod(l); p, type(p)
            (24, <class 'sage.interfaces.regina.ReginaElement'>)
            sage: c = CubicBraidGroup(3).an_element()
            sage: cr = regina(c); cr
            <regina.GroupExpression: g0 g1>
            sage: cr**(-3)
            <regina.GroupExpression: g1^-1 g0^-1 g1^-1 g0^-1 g1^-1 g0^-1>
        """
        P = self._check_valid()
        sinst = self._inst
        oinst = other if other is None else other._inst

        def is_native(inst):
            return type(inst) in (int, float, complex)
        if operation in ('+', '*'):
            if is_native(sinst) and is_native(oinst):
                if operation == '*':
                    return P(sinst * oinst)
                return P(sinst + oinst)
            if type(sinst) is type(oinst):
                if hasattr(self, 'addTermsLast'):
                    new = self.__deepcopy__()
                    new.addTermsLast(other)
                    new.simplify()
                    return new
        if operation == '^':
            if is_native(sinst) and is_native(oinst):
                return P(sinst**oinst)
            try:
                exp = int(other)
            except TypeError:
                raise TypeError('only integer exponents allowed!')

            if exp == 1:
                return self
            if exp == 2:
                return self * self
            if exp > 0:
                return self**(exp - 1) * self
            return (~self)**(-exp)
        if operation == '1/':
            if is_native(sinst):
                return P(1 / sinst)
            if hasattr(self, 'inverse'):
                return self.inverse()
        return super()._operation(operation, other=other)


@instancedoc
class PythonInternalFunctionElement(InterfaceFunctionElement):
    r"""
    Interface methods of interface elements.

    EXAMPLES::

        sage: A = regina.AbelianGroup()
        sage: type(A.addRank)
        <class 'sage.interfaces.python_internal.PythonInternalFunctionElement'>
        sage: M = snappy.Manifold('9_42')
        sage: M.DT_code
        <bound method Triangulation.DT_code of 9_42(0,0)>
        sage: type(M.DT_code)
        <class 'sage.interfaces.python_internal.PythonInternalFunctionElement'>
    """
    def __call__(self, *args, **kwds):
        r"""
        Call this function with the given args and kwds.

        This method is overloaded since the functions are
        Python functions which don't need to be interpreted
        from strings.

        EXAMPLES::

            sage: A = regina.AbelianGroup(); A
            <regina.AbelianGroup: 0>
            sage: a = A.addRank
            sage: a(); A
            <regina.AbelianGroup: Z>
            sage: a(int(2)); A
            <regina.AbelianGroup: 3 Z>
        """
        P = self._obj.parent()
        return P._function_call(self._name, *args, **kwds)


@instancedoc
class PythonInternalFunction(InterfaceFunction):
    r"""
    Interface Function.

    EXAMPLES::

        sage: m = regina.MatrixInt; m
        <class 'regina.engine.MatrixInt'>
        sage: type(m)
        <class 'sage.interfaces.python_internal.PythonInternalFunction'>
    """
    def __call__(self, *args, **kwds):
        r"""
        Call this function with the given args and kwds.

        This method is overloaded since the functions are
        Python functions which don't need to be interpreted
        from strings.

        EXAMPLES::

            sage: m = regina.MatrixInt
            sage: m([[1, 2], [3, 4]])
            <regina.MatrixInt: [[ 1 2 ] [ 3 4 ]]>
        """
        P = self._parent
        return P._function_call(self._name, *args, **kwds)
