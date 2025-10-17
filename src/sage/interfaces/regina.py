# sage.doctest: optional regina
r"""
Interface to Regina

`Regina <https://regina-normal.github.io/>`__ is an open source software for
low-dimensional topology. From the home-page:

.. NOTE::

    Regina is a software package for 3-manifold and 4-manifold
    topologists, with a focus on triangulations, knots and links,
    normal surfaces, and angle structures.

    For 3-manifolds, it includes high-level tasks such as 3-sphere
    recognition, connected sum decomposition and Hakenness testing,
    comes with a rich database of census manifolds, and incorporates
    the SnapPea kernel for working with hyperbolic manifolds.
    For 4-manifolds, it offers a range of combinatorial and algebraic
    tools, plus support for normal hypersurfaces. For knots and links,
    Regina can perform combinatorial manipulation, compute knot
    polynomials, and work with several import/export formats.

The Regina interface will only work if the optional Sage package Regina
is installed. The interface lets you send certain Sage objects to Regina,
run Regina functions, import certain Regina expressions to Sage,
or any combination of the above.

To send a Sage object ``sobj`` to Regina, call ``regina(sobj)``.
This exports the Sage object to Regina and returns a new Sage object
wrapping the Regina expression/variable, so that you can use the
Regina variable from within Sage. You can then call Regina
functions on the new object; for example::

    sage: F3 = FreeGroup(3)
    sage: F3r = regina(F3); F3r
    <regina.GroupPresentation: < a b c >>
    sage: (F3r.parent(), type(F3r))
    (Regina, <class 'sage.interfaces.regina.ReginaElement'>)
    sage: f = F3.an_element(); f
    x0*x1*x2
    sage: fr = regina(f); fr
    <regina.GroupExpression: g0 g1 g2>
    sage: rel = fr^2; rel
    <regina.GroupExpression: g0 g1 g2 g0 g1 g2>
    sage: (type(fr), type(rel))
    (<class 'sage.interfaces.regina.ReginaElement'>,
     <class 'sage.interfaces.regina.ReginaElement'>)
    sage: F3r.addRelation(rel); F3r
    <regina.GroupPresentation: < a b c | a b c a b c >>

In the above example the relations are added using Regina's
``addRelations`` method.

To see Regina's output you can simply print the Regina wrapper
object. However if you want to import Regina's output back to Sage,
call the Regina wrapper object's ``sage()`` method. This method returns
a native Sage object::

    sage: G3 = F3r.sage(); G3
    Finitely presented group < x0, x1, x2 | (x0*x1*x2)^2 >
    sage: type(G3)
    <class 'sage.groups.finitely_presented.FinitelyPresentedGroup_with_category'>
    sage: regina(G3) == F3r
    True
    sage: fr.sage() == f
    True
    sage: rel.sage() == f^2
    True

If you want to run a Regina function and don't already have the input
in the form of a Sage object, then it might be simpler to input a string
``expr`` to ``regina(expr)``. This string will be evaluated as if you had
typed it into Regina::

    sage: rL = regina("Link('dabcabcv-')"); rL
    <regina.Link: 3-crossing knot: +++ ( ^0 _1 ^2 _0 ^1 _2 )>

Alternatively, all constructors of Regina classes can be used directly as attributes
of the interface::

    sage: rL == regina.Link('dabcabcv-')
    True

Finally, if you just want to use a Regina command line from within
Sage, the IPython magic function ``%regina`` dumps you into an interactive
command-line Regina session. As long as you work in this environment the
prompt is ``regina:``. To finish the environment type ``CTRL+D``::

    sage: %regina                                 # not tested

    --> Switching to Regina <--

    regina: u = Link()
    None
    regina: u
    <regina.Link: Empty link>
    regina: type(u)
    <class 'regina.engine.Link'>
    regina: v = u.fromKnotSig('iabcdbefcdghaefghRsgF+m')
    None
    regina: v
    <regina.Link: 8-crossing knot: ++++--++ ( ^0 ^1 _2 ^3 _1 ^4 _5 ^2 _3 _6 ^7 _0 _4 ^5 ^6 _7 )>
    regina: type(v)
    <class 'regina.engine.Link'>
    regina: v.homfly()
    <regina.Laurent2: 2 x^-2 y^2 + 3 x^-2 - x^-4 y^4 - 3 x^-4 y^2 - 3 x^-4 + x^-6 y^2 + x^-6>
    regina: u.homfly()
    <regina.Laurent2: 0>
    regina: type(u)
    <class 'regina.engine.Link'>

    --> Exiting back to Sage <--

    sage:                                         # not tested


Complicated translations
------------------------

The ``robj.sage()`` method tries to convert a Regina object to a Sage
object. In many cases, it will just work. In particular, it should be able to
convert expressions entirely consisting of:

- numbers, i.e. integers, floats, complex numbers;
- functions and named constants also present in Sage, where:

    - Sage knows how to translate the function or constant's name from
      Regina's, or
    - the Sage name for the function or constant is trivially related to
      Regina's;

- symbolic variables whose names don't pathologically overlap with
  objects already defined in Sage.

This method will not work when Regina's output includes:

- strings;
- functions unknown to Sage;
- Regina functions with different parameters/parameter order to
  the Sage equivalent.


AUTHORS:

- Sebastian Oehms (2025): first version.
"""

##############################################################################
#       Copyright (C) 2025 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
##############################################################################

from enum import Enum
from sage.interfaces.interface import Interface, InterfaceElement, InterfaceFunction, InterfaceFunctionElement
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.misc.instancedoc import instancedoc


class AlgorithmExt(Enum):
    r"""
    Enum to select algorithm choices.

    This extends the ``Algorithm`` class of Regina.
    """
    ALG_SIMPLIFY = 4
    ALG_WIRTINGER = 5
    ALG_USE_EXTERIOR = 6


class Regina(ExtraTabCompletion, Interface):
    r"""
    Interface to the Regina interpreter.

    EXAMPLES::

        sage: K = Knots().from_table(8, 21)
        sage: Kr = regina(K); Kr
        <regina.Link: 8-crossing knot: ----++-- ( _5 _0 ^1 _2 _3 ^6 _7 ^3 _4 ^5 _6 ^7 ^0 _1 ^2 ^4 )>
        sage: Kr.knotSig()
        'iabcdbefcdghaefghRsgF+m'

    More examples can be found in the module header.
    """
    def __init__(self):
        r"""
        Python constructor.

        TESTS::

            sage: TestSuite(regina).run(skip=['_test_pickling', '_test_category'])
        """
        Interface.__init__(self, name='regina')
        self._initialized = False  # done lazily
        self._regina_globals = None
        self._namespace = None

    def _lazy_init(self):
        r"""
        Initialize the Regina interpreter.

        Implemented according to R interface.

        EXAMPLES::

            sage: regina._lazy_init()
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
        Start up the Regina interpreter and sets the initial prompt and options.

        This is called the first time the Regina interface is actually used.

        EXAMPLES::

            sage: regina._start()
            sage: type(regina._namespace.Cyclotomic())
            <class 'regina.engine.Cyclotomic'>
        """
        if not self._regina_globals:
            from sage.features.interfaces import Regina
            Regina().module.require()
            import regina
            self._namespace = regina.engine
            d = self._namespace.__dict__
            # add extras to the fixed namespace
            d['AlgorithmExt'] = AlgorithmExt
            for e in AlgorithmExt:
                d[e.name] = e
            # set up the globals
            D = {k: v for k, v in d.items()}
            D[self._namespace.__name__] = self._namespace
            D[self.name()] = regina
            self._regina_globals = D

    def _install_hints(self):
        """
        Hints for installing regina on your computer.

        EXAMPLES::

            sage: len(regina._install_hints())
            99
        """
        return """
In order to use the Regina interface you need to have the
optional Sage package Regina installed.
"""

    def _eval(self, code):
        """
        Evaluates a command inside the Regina interpreter and returns the output
        as a Regina result.

        EXAMPLES::

            sage: regina._eval("Link('iabcdbefcdghaefghRsgF+m')")
            <regina.Link: 8-crossing knot: ++++--++ ( ^0 ^1 _2 ^3 _1 ^4 _5 ^2 _3 _6 ^7 _0 _4 ^5 ^6 _7 )>
        """
        self._lazy_init()
        globs = self._regina_globals
        if type(code) is str and code.find('=') < 0:
            pre_bracket = code.split('(')[0]
            if pre_bracket in globs:
                val = globs[pre_bracket]
                nam = self._namespace.__name__
                if callable(val):
                    return eval('%s.%s' % (nam, code), globs)
                else:
                    return val
            else:
                return eval(code, globs)
        return exec(code, globs)

    def eval(self, code, *args, **kwds):
        """
        Evaluates a command inside the Regina interpreter and returns the output
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
        """
        return self._regina_globals[var]

    def set(self, var, value):
        """
        Set the variable ``var`` to the given ``value``.

        EXAMPLES::

            sage: regina.set('myLink', 'Link')
            sage: regina.get('myLink')
            <class 'regina.engine.Link'>
        """
        self._lazy_init()
        globs = self._regina_globals
        if not isinstance(value, str):
            globs[var] = value
            return
        try:
            val = self._eval(value)
            globs[var] = val
        except (NameError, AttributeError, KeyError):
            pass
        super().set(var, value)

    def _regina_object(self, x) -> bool:
        r"""
        Return ``True`` if ``x`` is an instance of a class
        from Regina's namespace.

        EXAMPLES::

            sage: L = regina.get('Link')
            sage: regina._regina_object(L)
            False
            sage: regina._regina_object(L())
            True
            sage: regina._regina_object([])
            False
        """
        self._lazy_init()
        d = self._namespace.__dict__
        cl = x.__class__
        cln = cl.__name__
        if cln in d:
            if cl == d[cln]:
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

        This method is overloaded to add Python types from regina._namespace

        EXAMPLES::

            sage: L = regina.get('Link')
            sage: regina._coerce_impl(L())
            <regina.Link: Empty link>
        """
        if self._regina_object(x):
            return self(self._create(x))
        return super()._coerce_impl(x, use_special=use_special)

    def _convert_args_kwds(self, *args, **kwds):
        """
        Convert all of the ``args`` and ``kwds`` to instances of Regina
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
            if isinstance(arg, InterfaceElement) and arg.parent() is self:
                return arg._inst
            elif isinstance(arg, (list, tuple)):
                return type(arg)([convert_arg(i) for i in arg])
            elif hasattr(arg, '_regina_'):
                reg = arg._regina_(self)
                return convert_arg(reg)
            else:
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
            <regina.Polynomial: 5/3 x - 3>
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
        else:
            if len(kwds) == 0:
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
            if self._regina_object(res):
                new = res.__class__(res)  # this is the way to get a copy of a Regina object
                return self(new)
            else:
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
        """
        return ReginaElement

    def help(self, cmd, long=False):
        r"""
        Return the Regina documentation of the given command.

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

           Currently returns all keys of the namspace dictionary.

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
            <class 'sage.interfaces.regina.ReginaFunction'>
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
            return ReginaFunction(self, self._namespace.__dict__[attrname])
        res = self(self._create(attr))
        res._inst = attr
        return res


@instancedoc
class ReginaElement(ExtraTabCompletion, InterfaceElement):
    r"""
    Element class of the Regina interface.

    Its instances are usually constructed via the instance call of its parent.
    It wrapes the Regina library for this object. In a session Regina methods
    can be obtained using tab completion.

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

    Access to the Regina expression objects::

        sage: res = re._inst
        sage: type(res)
         <class 'regina.engine.GroupExpression'>

    Applying Regina methods::

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
        return self.parent().new(self._inst[n])

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
            <class 'sage.interfaces.regina.ReginaFunctionElement'>
            sage: regina.AbelianGroup().detail._name
            <bound method pybind11_detail_function_record_v1_system_libstdcpp_gxx_abi_1xxx_use_cxx11_abi_0.detail of <regina.AbelianGroup: 0>>
        """
        P = self._check_valid()
        if attrname == '_inst':
            self._inst = P.get(self.name())
            return self._inst
        elif attrname[:1] == "_":
            raise AttributeError
        else:
            inst = self._inst
            if hasattr(inst, attrname):
                attr = inst.__getattribute__(attrname)
                if callable(attr):
                    return ReginaFunctionElement(self, attr)
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

    def __deepcopy__(self, memo=None):
        r"""
        EXAMPLES::

            sage: C = CubicBraidGroup(3)
            sage: c = C.an_element()
            sage: cr = regina(c); cr
            <regina.GroupExpression: g0 g1>
            sage: from copy import deepcopy, copy
            sage: crd = deepcopy(cr)
            sage: crc = copy(cr)
            sage: cr.cycleRight(); cr
            <regina.GroupExpression: g1 g0>
            sage: (crd, crd == cr)
            (<regina.GroupExpression: g0 g1>, False)
            sage: (crc, crc == cr)
            (<regina.GroupExpression: g1 g0>, True)
        """
        P = self._check_valid()
        inst = self._inst
        new = inst.__class__(inst)
        res = P(new)
        res._sage_parent = self._sage_parent
        return res

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
        from sage.structure.richcmp import rich_to_bool, op_EQ, op_NE
        if self._inst == other._inst:
            return rich_to_bool(op, 0)
        elif op == op_EQ:
            return False
        elif op == op_NE:
            return True
        try:
            if self._inst < other._inst:
                return rich_to_bool(op, -1)
            elif self._inst > other._inst:
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
                else:
                    return P(sinst + oinst)
            if type(sinst) == type(oinst):
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
            elif exp == 2:
                return self * self
            elif exp > 0:
                for i in range(exp):
                    return self**(exp - 1) * self
            else:
                return (~self)**(-exp)
        if operation == '1/':
            if is_native(sinst):
                return P(1 / sinst)
            if hasattr(self, 'inverse'):
                return self.inverse()
        return super()._operation(operation, other=other)

    def _sage_(self, locals={}):
        r"""
        Attempt to return a Sage version of this object.

        This method works successfully when Regina returns a result
        or list of results that consist only of:

        - numbers, i.e. integers, floats, complex numbers;
        - functions and named constants also present in Sage, where:

          * Sage knows how to translate the function or constant's name
            from Regina's naming scheme, or
          * you provide a translation dictionary `locals`, or
          * the Sage name for the function or constant is simply the
            Regina name in lower case;

        - symbolic variables whose names do not pathologically overlap with
          objects already defined in Sage.

        This method will not work when Regina's output includes:

        - strings;
        - functions unknown to Sage;
        - Regina functions with different parameters/parameter order to
          the Sage equivalent. In this case, define a function to do the
          parameter conversion, and pass it in via the locals dictionary.

        EXAMPLES::

            sage: p = regina("Laurent2([(1,2,3), (2,-1,5)])"); p
            <regina.Laurent2: 5 x^2 y^-1 + 3 x y^2>
            sage: p.sage()
            3*x*y^2 + 5*x^2*y^-1
            sage: R.<u> = PolynomialRing(ZZ)
            sage: p = u**3 -2*u + 9
            sage: rp = regina(p); rp
            <regina.Polynomial: x^3 - 2 x + 9>
            sage: rp.sage() == p
            True
            sage: F3 = FreeGroup(3)
            sage: a, b, c = F3.gens()
            sage: f = a**2*b*~c
            sage: fr = regina(f); fr
            <regina.GroupExpression: g0^2 g1 g2^-1>
            sage: fr.sage() == f
            True
        """
        def from_detail_str(lc):
            r"""
            Regina provides a detail method for many of its classes.
            Here we try to use it to convert back to Sage.
            """
            if locals:
                lc.update(locals)
            from sage.misc.sage_eval import sage_eval
            s = self.detail().split('\n')[0]
            s = s.replace(' ', '')
            v = list(lc)
            for i in v + ['(']:
                for j in v + [')']:
                    s = s.replace('%s%s' % (i, j), '%s*%s' % (i, j))
            s = implicit_mul(s)
            return sage_eval(s, locals=lc)

        P = self._check_valid()
        inst = self._inst
        nspc = P._namespace
        if hasattr(inst, 'sage'):
            return inst.sage()
        elif isinstance(inst, (nspc.Polynomial, nspc.Laurent, nspc.Laurent2)):
            if self._sage_parent:
                R = self._sage_parent
                old_var_names = ['x', 'y']
                new_var_names = R.variable_names()
                lc = {old_var_names[i]: R.gens_dict()[new_var_names[i]] for i in range(len(new_var_names))}
            elif isinstance(inst, nspc.Polynomial):
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                from sage.rings.integer_ring import ZZ
                R = PolynomialRing(ZZ, 'x')
                lc = R.gens_dict()
            else:
                from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
                from sage.rings.integer_ring import ZZ
                if isinstance(inst, nspc.Laurent):
                    R = LaurentPolynomialRing(ZZ, 'x')
                else:
                    R = LaurentPolynomialRing(ZZ, 'x, y')
                lc = R.gens_dict()
            return from_detail_str(lc)
        elif isinstance(inst, nspc.GroupExpression):
            from sage.repl.preparse import implicit_mul
            num_gens = max(t.generator for t in inst.terms()) + 1
            if self._sage_parent:
                F = self._sage_parent
            else:
                from sage.groups.free_group import FreeGroup
                F = FreeGroup(num_gens)
            gens = F.gens()
            lc = {'g%s' % i: gens[i] for i in range(num_gens)}
            return from_detail_str(lc)
        elif isinstance(inst, nspc.Link):
            from sage.knots.link import Link
            return Link(inst.pdData())
        elif hasattr(self, 'detail'):
            return from_detail_str(locals)
        elif locals:
            # if locals are given we use `_sage_repr`
            # surely this only covers simple cases
            from sage.misc.sage_eval import sage_eval
            return sage_eval(self._sage_repr(), locals=locals)
        return inst


@instancedoc
class ReginaFunctionElement(InterfaceFunctionElement):
    r"""
    Interface methods of interface elements.

    EXAMPLES::

        sage: A = regina.AbelianGroup()
        sage: A.addRank
        <bound method pybind11_detail_function_record_v1_system_libstdcpp_gxx_abi_1xxx_use_cxx11_abi_0.addRank of <regina.AbelianGroup: 0>>
        sage: type(A.addRank)
        <class 'sage.interfaces.regina.ReginaFunctionElement'>
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
class ReginaFunction(InterfaceFunction):
    r"""
    Interface Function.

    EXAMPLES::

        sage: m = regina.MatrixInt; m
        <class 'regina.engine.MatrixInt'>
        sage: type(m)
        <class 'sage.interfaces.regina.ReginaFunction'>
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


# An instance
regina = Regina()
