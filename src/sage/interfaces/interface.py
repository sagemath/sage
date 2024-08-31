# sage.doctest: needs sage.libs.gap sage.libs.pari sage.libs.singular sage.symbolic
r"""
Common Interface Functionality

See the examples in the other sections for how to use specific
interfaces. The interface classes all derive from the generic
interface that is described in this section.

AUTHORS:

- William Stein (2005): initial version

- William Stein (2006-03-01): got rid of infinite loop on startup if
  client system missing

- Felix Lawrence (2009-08-21): edited ._sage_() to support lists and float exponents in foreign notation.

- Simon King (2010-09-25): Expect._local_tmpfile() depends on
  Expect.pid() and is cached; Expect.quit() clears that cache,
  which is important for forking.

- Jean-Pierre Flori (2010,2011): Split non Pexpect stuff into a parent class.

- Simon King (2015): Improve pickling for InterfaceElement
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator

from sage.structure.sage_object import SageObject
from sage.structure.parent_base import ParentWithBase
from sage.structure.element import Element, parent
from sage.structure.richcmp import rich_to_bool

import sage.misc.sage_eval
from sage.misc.fast_methods import WithEqualityById
from sage.misc.instancedoc import instancedoc


class AsciiArtString(str):
    def __repr__(self):
        return str(self)


class Interface(WithEqualityById, ParentWithBase):
    """
    Interface interface object.

    .. NOTE::

        Two interfaces compare equal if and only if they are identical
        objects (this is a critical constraint so that caching of
        representations of objects in interfaces works
        correctly). Otherwise they are never equal.
    """
    def __init__(self, name):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Maxima() == maxima
            False
            sage: maxima == maxima
            True

            sage: Maxima() != maxima
            True
            sage: maxima != maxima
            False
        """
        self.__name = name
        self.__coerce_name = '_' + name.lower() + '_'
        self.__seq = -1
        self._available_vars = []
        self._seed = None
        ParentWithBase.__init__(self, self)

    def _repr_(self):
        return self.__name.capitalize()

    def name(self, new_name=None):
        return self.__name

    def get_seed(self):
        """
        Return the seed used to set the random number generator in
        this interface.

        The seed is initialized as ``None`` but should be set when the
        interface starts.

        EXAMPLES::

            sage: s = Singular()
            sage: s.set_seed(107)
            107
            sage: s.get_seed()
            107
        """
        return self._seed

    def rand_seed(self):
        """
        Return a random seed that can be put into ``set_seed`` function
        for any interpreter.

        This should be overridden if the particular interface needs
        something other than a small positive integer.

        EXAMPLES::

            sage: from sage.interfaces.interface import Interface
            sage: i = Interface("")
            sage: i.rand_seed() # random
            318491487

            sage: s = Singular()
            sage: s.rand_seed() # random
            365260051
        """
        import sage.doctest
        if sage.doctest.DOCTEST_MODE:
            # set the random seed through the current randstate
            from sage.misc.randstate import current_randstate
            seed = current_randstate().seed()
        else:
            from sage.misc.randstate import randstate
            seed = randstate().seed()

        return seed & 0x1FFFFFFF

    def set_seed(self, seed=None):
        """
        Set the random seed for the interpreter and return the new
        value of the seed.

        This is dependent on which interpreter so must be implemented
        in each separately. For examples see gap.py or singular.py.

        If seed is ``None`` then should generate a random seed.

        EXAMPLES::

            sage: s = Singular()
            sage: s.set_seed(1)
            1
            sage: [s.random(1,10) for i in range(5)]
            [8, 10, 4, 9, 1]

            sage: from sage.interfaces.interface import Interface
            sage: i = Interface("")
            sage: i.set_seed()
            Traceback (most recent call last):
            ...
            NotImplementedError: This interpreter did not implement a set_seed function
        """
        raise NotImplementedError("This interpreter did not implement a set_seed function")

    def interact(self):
        r"""
        This allows you to interactively interact with the child
        interpreter.

        Press :kbd:`Ctrl` + :kbd:`D` or type 'quit' or 'exit' to exit and
        return to Sage.

        .. NOTE::

           This is completely different than the console() member
           function. The console function opens a new copy of the
           child interpreter, whereas the interact function gives you
           interactive access to the interpreter that is being used by
           Sage. Use sage(xxx) or interpretername(xxx) to pull objects
           in from sage to the interpreter.
        """
        from sage.repl.interpreter import interface_shell_embed
        shell = interface_shell_embed(self)
        try:
            ipython = get_ipython()
        except NameError:
            shell()
        else:
            shell(local_ns=dict(ipython.user_ns))

    def _pre_interact(self):
        pass

    def _post_interact(self):
        pass

    def cputime(self):
        """
        CPU time since this process started running.
        """
        raise NotImplementedError

    def read(self, filename):
        r"""
        EXAMPLES::

            sage: filename = tmp_filename()
            sage: f = open(filename, 'w')
            sage: _ = f.write('x = 2\n')
            sage: f.close()
            sage: octave.read(filename)  # optional - octave
            sage: octave.get('x')        # optional - octave
            ' 2'
            sage: import os
            sage: os.unlink(filename)
        """
        self.eval(self._read_in_file_command(filename))

    def _read_in_file_command(self, filename):
        raise NotImplementedError

    def eval(self, code, **kwds):
        """
        Evaluate code in an interface.

        This method needs to be implemented in sub-classes.

        Note that it is not always to be expected that
        it returns a non-empty string. In contrast,
        :meth:`get` is supposed to return the result of applying
        a print command to the object so that the output is easier
        to parse.

        Likewise, the method :meth:`_eval_line` for evaluation of a single
        line, often makes sense to be overridden.
        """
        raise NotImplementedError

    _eval_line = eval

    def execute(self, *args, **kwds):
        return self.eval(*args, **kwds)

    def __call__(self, x, name=None):
        r"""
        Create a new object in ``self`` from ``x``.

        The object ``X`` returned can be used like any Sage object, and
        wraps an object in ``self``.  The standard arithmetic operators
        work.  Moreover if ``foo`` is a function then::

            ``X.foo(y,z,...)``

        calls ``foo(X, y, z, ...)`` and returns the corresponding object.

        EXAMPLES::

            sage: gp(2)
            2
            sage: gp('2')
            2
            sage: a = gp(2); gp(a) is a
            True

        TESTS:

        Check conversion of Booleans (:issue:`28705`)::

            sage: giac(True)
            true
            sage: maxima(True)
            true
        """
        cls = self._object_class()

        # Handle the case when x is an object
        # in some interface.
        if isinstance(x, InterfaceElement):
            if x.parent() is self:
                return x

            # We convert x into an object in this
            # interface by first going through Sage.
            try:
                return self(x._sage_())
            except (NotImplementedError, TypeError):
                pass

        if isinstance(x, str):
            return cls(self, x, name=name)
        try:
            # Special methods do not and should not have an option to
            # set the name directly, as the identifier assigned by the
            # interface should stay consistent. An identifier with a
            # user-assigned name might change its value, so we return a
            # new element.
            result = self._coerce_from_special_method(x)
            return result if name is None else result.name(new_name=name)
        except TypeError:
            raise
        except AttributeError:
            pass
        try:
            result = self._coerce_impl(x, use_special=False)
            return result if name is None else result.name(new_name=name)
        except TypeError as msg:
            try:
                return cls(self, str(x), name=name)
            except TypeError:
                raise TypeError(msg)

    def _coerce_from_special_method(self, x):
        """
        Try to coerce to ``self`` by calling a special underscore method.

        If no such method is defined, raises an :exc:`AttributeError` instead
        of a :exc:`TypeError`.
        """
        s = '_%s_' % self.name()
        if s == '_maxima_lib_':
            s = '_maxima_'
        if s == '_pari_':
            s = '_gp_'
        try:
            return (x.__getattribute__(s))(self)
        except AttributeError:
            return self(x._interface_init_())

    def _coerce_impl(self, x, use_special=True):
        r"""
        Coerce pure Python types via corresponding Sage objects.

        TESTS:

        Check that python type ``complex`` can be converted (:issue:`31775`)::

            sage: giac(complex(I))**2  # should not return `j^2`
            -1
        """
        if isinstance(x, bool):
            return self(self._true_symbol() if x else self._false_symbol())
        elif isinstance(x, int):
            from sage.rings.integer import Integer
            return self(Integer(x))
        elif isinstance(x, float):
            from sage.rings.real_double import RDF
            return self(RDF(x))
        elif isinstance(x, complex):
            from sage.rings.complex_double import CDF
            return self(CDF(x))
        if use_special:
            try:
                return self._coerce_from_special_method(x)
            except AttributeError:
                pass

        if isinstance(x, (list, tuple)):
            A = []
            z = []
            cls = self._object_class()
            for v in x:
                if isinstance(v, cls):
                    A.append(v.name())
                    z.append(v)
                else:
                    w = self(v)
                    A.append(w.name())
                    z.append(w)
            X = ','.join(A)
            r = self.new('%s%s%s' % (self._left_list_delim(), X, self._right_list_delim()))
            r.__sage_list = z   # do this to avoid having the entries of the list be garbage collected
            return r

        raise TypeError("unable to coerce element into %s" % self.name())

    def new(self, code):
        return self(code)

    ###################################################################
    # these should all be appropriately overloaded by the derived class
    ###################################################################

    def _left_list_delim(self):
        return "["

    def _right_list_delim(self):
        return "]"

    def _left_func_delim(self):
        return "("

    def _right_func_delim(self):
        return ")"

    def _assign_symbol(self):
        return "="

    def _equality_symbol(self):
        raise NotImplementedError

    # For efficiency purposes, you should definitely override these
    # in your derived class.
    def _true_symbol(self):
        try:
            return self.__true_symbol
        except AttributeError:
            self.__true_symbol = self.get('1 %s 1' % self._equality_symbol())
            return self.__true_symbol

    def _false_symbol(self):
        try:
            return self.__false_symbol
        except AttributeError:
            self.__false_symbol = self.get('1 %s 2' % self._equality_symbol())
            return self.__false_symbol

    def _lessthan_symbol(self):
        return '<'

    def _greaterthan_symbol(self):
        return '>'

    def _inequality_symbol(self):
        return '!='

    def _relation_symbols(self):
        """
        Return a dictionary with operators as the keys and their
        string representation as the values.

        EXAMPLES::

            sage: import operator
            sage: symbols = mathematica._relation_symbols()
            sage: symbols[operator.eq]
            '=='
        """
        return dict([(operator.eq, self._equality_symbol()), (operator.ne, self._inequality_symbol()),
                     (operator.lt, self._lessthan_symbol()), (operator.le, "<="),
                     (operator.gt, self._greaterthan_symbol()), (operator.ge, ">=")])

    def _exponent_symbol(self):
        """
        Return the symbol used to denote ``*10^`` in floats, e.g 'e' in 1.5e6.

        EXAMPLES::

            sage: from sage.interfaces.expect import Expect
            sage: Expect('nonexistent_interface', 'fake')._exponent_symbol()
            'e'
        """
        return 'e'

    ############################################################
    #         Functions for working with variables.
    #  The first three must be overloaded by derived classes,
    #  and the definition depends a lot on the class.  But
    #  the functionality one gets from this is very nice.
    ############################################################

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '%s%s%s;' % (var, self._assign_symbol(), value)
        self.eval(cmd)

    def get(self, var):
        """
        Get the value of the variable var.

        Note that this needs to be overridden in some interfaces,
        namely when getting the string representation of an object
        requires an explicit print command.
        """
        return self.eval(var)

    def get_using_file(self, var):
        r"""
        Return the string representation of the variable var in self,
        possibly using a file. Use this if var has a huge string
        representation, since it may be way faster.

        .. warning::

           In fact unless a special derived class implements this, it
           will *not* be any faster. This is the case for this class
           if you're reading it through introspection and seeing this.
        """
        return self.get(var)

    def clear(self, var):
        """
        Clear the variable named var.
        """
        self._available_vars.append(var)

    def _default_var_name(self):
        return "sage"

    def _next_var_name(self):
        if len(self._available_vars) != 0:
            v = self._available_vars[0]
            del self._available_vars[0]
            return v
        self.__seq += 1
        return self._default_var_name() + str(self.__seq)

    def _create(self, value, name=None):
        name = self._next_var_name() if name is None else name
        self.set(name, value)
        return name

    def _object_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.expect import Expect
            sage: Expect._object_class(maxima)
            <class 'sage.interfaces.expect.ExpectElement'>
        """
        return InterfaceElement

    def _function_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.interface import Interface
            sage: Interface._function_class(maxima)
            <class 'sage.interfaces.interface.InterfaceFunction'>
        """
        return InterfaceFunction

    def _function_element_class(self):
        """
        EXAMPLES::

            sage: from sage.interfaces.interface import Interface
            sage: Interface._function_element_class(maxima)
            <class 'sage.interfaces.interface.InterfaceFunctionElement'>
        """
        return InterfaceFunctionElement

    def _convert_args_kwds(self, args=None, kwds=None):
        """
        Convert all of the ``args`` and ``kwds`` to be elements of this
        interface.

        EXAMPLES::

            sage: args = [5]
            sage: kwds = {'x': 6}
            sage: args, kwds = gap._convert_args_kwds(args, kwds)
            sage: args
            [5]
            sage: list(map(type, args))
            [<class 'sage.interfaces.gap.GapElement'>]
            sage: type(kwds['x'])
            <class 'sage.interfaces.gap.GapElement'>
        """
        args = [] if args is None else args
        kwds = {} if kwds is None else kwds
        if not isinstance(args, list):
            args = [args]
        for i, arg in enumerate(args):
            if not isinstance(arg, InterfaceElement) or arg.parent() is not self:
                args[i] = self(arg)
        for key, value in kwds.items():
            if not isinstance(value, InterfaceElement) or value.parent() is not self:
                kwds[key] = self(value)

        return args, kwds

    def _check_valid_function_name(self, function):
        """
        Check to see if function is a valid function name in this
        interface. If it is not, an exception is raised. Otherwise, nothing
        is done.

        EXAMPLES::

            sage: gap._check_valid_function_name('SymmetricGroup')
            sage: gap._check_valid_function_name('')
            Traceback (most recent call last):
            ...
            ValueError: function name must be nonempty
            sage: gap._check_valid_function_name('__foo')
            Traceback (most recent call last):
            ...
            AttributeError
        """
        if function == '':
            raise ValueError("function name must be nonempty")
        if function[:2] == "__":
            raise AttributeError

    def function_call(self, function, args=None, kwds=None):
        """
        EXAMPLES::

            sage: maxima.quad_qags(x, x, 0, 1, epsrel=1e-4)
            [0.5,5.5511151231257...e-15,21,0]
            sage: maxima.function_call('quad_qags', [x, x, 0, 1], {'epsrel':'1e-4'})
            [0.5,5.5511151231257...e-15,21,0]
        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)
        s = self._function_call_string(function,
                                       [s.name() for s in args],
                                       ['%s=%s' % (key, value.name()) for key, value in kwds.items()])
        return self.new(s)

    def _function_call_string(self, function, args, kwds):
        """
        Return the string used to make function calls.

        EXAMPLES::

            sage: maxima._function_call_string('diff', ['f(x)', 'x'], [])
            'diff(f(x),x)'
        """
        return "%s(%s)" % (function, ",".join(list(args) + list(kwds)))

    def call(self, function_name, *args, **kwds):
        return self.function_call(function_name, args, kwds)

    def _contains(self, v1, v2):
        raise NotImplementedError

    def __getattr__(self, attrname):
        """
        TESTS::

            sage: from sage.structure.parent_base import ParentWithBase
            sage: from sage.interfaces.singular import singular
            sage: ParentWithBase.__getattribute__(singular, '_coerce_map_from_')
            <bound method Singular._coerce_map_from_ of Singular>
        """
        try:
            return ParentWithBase.__getattribute__(self, attrname)
        except AttributeError:
            if attrname[:1] == "_":
                raise
            return self._function_class()(self, attrname)

    def console(self):
        raise NotImplementedError

    def help(self, s):
        return AsciiArtString('No help on %s available' % s)


@instancedoc
class InterfaceFunction(SageObject):
    """
    Interface function.
    """
    def __init__(self, parent, name):
        self._parent = parent
        self._name = name

    def _repr_(self):
        return "%s" % self._name

    def __call__(self, *args, **kwds):
        return self._parent.function_call(self._name, list(args), kwds)

    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: gp.gcd.__doc__
            'gcd(x,{y}): greatest common divisor of x and y.'
        """
        M = self._parent
        return M.help(self._name)


@instancedoc
class InterfaceFunctionElement(SageObject):
    """
    Interface function element.
    """
    def __init__(self, obj, name):
        self._obj = obj
        self._name = name

    def _repr_(self):
        return "%s" % self._name

    def __call__(self, *args, **kwds):
        return self._obj.parent().function_call(self._name, [self._obj] + list(args), kwds)

    def help(self):
        print(self.__doc__)

    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: gp(2).gcd.__doc__
            'gcd(x,{y}): greatest common divisor of x and y.'
        """
        M = self._obj.parent()
        return M.help(self._name)


def is_InterfaceElement(x):
    """
    Return ``True`` if ``x`` is of type :class:`InterfaceElement`.

    EXAMPLES::

        sage: from sage.interfaces.interface import is_InterfaceElement
        sage: is_InterfaceElement(2)
        doctest:...: DeprecationWarning: the function is_InterfaceElement is deprecated; use isinstance(x, sage.interfaces.abc.InterfaceElement) instead
        See https://github.com/sagemath/sage/issues/34804 for details.
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(34804, "the function is_InterfaceElement is deprecated; use isinstance(x, sage.interfaces.abc.InterfaceElement) instead")

    return isinstance(x, InterfaceElement)


@instancedoc
class InterfaceElement(Element):
    """
    Interface element.
    """
    def __init__(self, parent, value, is_name=False, name=None):
        Element.__init__(self, parent)
        self._create = value
        if parent is None:
            return     # means "invalid element"
        # idea: Joe Wetherell -- try to find out if the output
        # is too long and if so get it using file, otherwise
        # don't.

        if is_name:
            self._name = value
        else:
            try:
                self._name = parent._create(value, name=name)
            except (TypeError, RuntimeError, ValueError) as x:
                raise TypeError(x)

    def _latex_(self):
        #        return "\\begin{verbatim}%s\\end{verbatim}"%self
        string = str(self)
        if '|' not in string:
            delim = '|'
        elif '#' not in string:
            delim = '#'
        elif '@' not in string:
            delim = '@'
        elif '~' not in string:
            delim = '~'
        return "\\verb%s%s%s" % (delim, string, delim)

    def __iter__(self):
        for i in range(1, len(self) + 1):
            yield self[i]

    def __len__(self):
        """
        Call self.sage() and return the length of that sage object.

        This approach is inefficient - each interface should override
        this method with one that calls the external program's length
        function.

        EXAMPLES::

            sage: len(gp([1,2,3]))
            3

        AUTHORS:

        - Felix Lawrence (2009-08-21)
        """
        return len(self.sage())

    def __reduce__(self):
        """
        The default linearisation is to return ``self``'s parent,
        which will then get the items returned by :meth:`_reduce`
        as arguments to reconstruct the element.

        EXAMPLES::

            sage: G = gap.SymmetricGroup(6)
            sage: loads(dumps(G)) == G     # indirect doctest
            True
            sage: y = gap(34)
            sage: loads(dumps(y))
            34
            sage: type(_)
            <class 'sage.interfaces.gap.GapElement'>
            sage: y = singular(34)
            sage: loads(dumps(y))
            34
            sage: type(_)
            <class 'sage.interfaces.singular.SingularElement'>
            sage: G = gap.PolynomialRing(QQ, ['x'])
            sage: loads(dumps(G))
            PolynomialRing( Rationals, ["x"] )
            sage: S = singular.ring(0, ('x'))
            sage: loads(dumps(S))
            polynomial ring, over a field, global ordering
            //   coefficients: QQ
            //   number of vars : 1
            //        block   1 : ordering lp
            //                  : names    x
            //        block   2 : ordering C

        Here are further examples of pickling of interface elements::

            sage: loads(dumps(gp('"abc"')))
            abc
            sage: loads(dumps(gp([1,2,3])))
            [1, 2, 3]
            sage: loads(dumps(pari('"abc"')))
            "abc"
            sage: loads(dumps(pari([1,2,3])))
            [1, 2, 3]
            sage: loads(dumps(r('"abc"')))                                        # optional - rpy2
            [1] "abc"
            sage: loads(dumps(r([1,2,3])))                                        # optional - rpy2
            [1] 1 2 3
            sage: loads(dumps(maxima([1,2,3])))
            [1,2,3]

        Unfortunately, strings in maxima can't be pickled yet::

            sage: loads(dumps(maxima('"abc"')))
            Traceback (most recent call last):
            ...
            TypeError: unable to make sense of Maxima expression '"abc"' in Sage
        """
        return self.parent(), (self._reduce(),)

    def _reduce(self):
        """
        Helper for pickling.

        By default, if ``self`` is a string, then the representation of
        that string is returned (not the string itself). Otherwise,
        it is attempted to return the corresponding Sage object.
        If this fails with a :exc:`NotImplementedError`, the string
        representation of ``self`` is returned instead.

        EXAMPLES::

            sage: S = singular.ring(0, ('x'))
            sage: S._reduce()
            Univariate Polynomial Ring in x over Rational Field
            sage: G = gap.PolynomialRing(QQ, ['x'])
            sage: G._reduce()
            'PolynomialRing( Rationals, ["x"] )'
            sage: G.sage()
            Traceback (most recent call last):
            ...
            NotImplementedError: Unable to parse output: PolynomialRing( Rationals, ["x"] )
            sage: singular('"abc"')._reduce()
            "'abc'"
            sage: singular('1')._reduce()
            1

        TESTS:

        Special care has to be taken with strings. Since for example `r("abc")` will be
        interpreted as the R-command abc (not a string in R), we have to reduce to
        `"'abc'"` instead. That is dependant on the Elements `is_string` function to
        be implemented correctly. This has gone wrong in the past and remained uncaught
        by the doctests because the original identifier was reused. This test makes sure
        that does not happen again::

            sage: # optional - rpy2
            sage: a = r("'abc'")
            sage: b = dumps(a)
            sage: r.set(a.name(), 0) # make sure that identifier reuse
            ....:                    # does not accidentally lead to success
            sage: loads(b)
            [1] "abc"
        """
        if self.is_string():
            return repr(self.sage())
        try:
            return self.sage()
        except NotImplementedError:
            return repr(self)

    def __call__(self, *args):
        self._check_valid()
        P = self.parent()
        return getattr(P, self.name())(*args)

    def __contains__(self, x):
        P = self._check_valid()
        if not isinstance(x, InterfaceElement) or x.parent() is not self.parent():
            x = P.new(x)
        return P._contains(x.name(), self.name())

    def _instancedoc_(self):
        """
        EXAMPLES::

            sage: gp(2).__doc__
            '2'
        """
        return str(self)

    def __hash__(self):
        """
        Return the hash of ``self``. This is a default implementation of hash
        which just takes the hash of the string of ``self``.
        """
        return hash('%s' % self)

    def _richcmp_(self, other, op):
        """
        Comparison of interface elements.

        NOTE:

        GAP has a special role here. It may in some cases raise an error
        when comparing objects, which is unwanted in Python. We catch
        these errors. Moreover, GAP does not recognise certain objects as
        equal even if there definitions are identical.

        NOTE:

        This methods need to be overridden if the subprocess would
        not return a string representation of a boolean value unless
        an explicit print command is used.

        TESTS:

        Here are examples in which GAP succeeds with a comparison::

            sage: gap('SymmetricGroup(8)')==gap('SymmetricGroup(8)')
            True
            sage: gap('SymmetricGroup(8)')>gap('AlternatingGroup(8)')
            False
            sage: gap('SymmetricGroup(8)')<gap('AlternatingGroup(8)')
            True

        Here, GAP fails to compare, and so ``False`` is returned.
        In previous Sage versions, this example actually resulted
        in an error; compare :issue:`5962`.
        ::

            sage: gap('DihedralGroup(8)')==gap('DihedralGroup(8)')
            False
        """
        P = self._check_valid()
        try:
            if P.eval("%s %s %s" % (self.name(), P._equality_symbol(),
                                     other.name())) == P._true_symbol():
                return rich_to_bool(op, 0)
        except RuntimeError:
            pass
        try:
            if P.eval("%s %s %s" % (self.name(), P._lessthan_symbol(), other.name())) == P._true_symbol():
                return rich_to_bool(op, -1)
        except RuntimeError:
            pass
        try:
            if P.eval("%s %s %s" % (self.name(), P._greaterthan_symbol(), other.name())) == P._true_symbol():
                return rich_to_bool(op, 1)
        except Exception:
            pass

        return NotImplemented

    def is_string(self):
        """
        Tell whether this element is a string.

        By default, the answer is negative.
        """
        return False

    def _matrix_(self, R):
        raise NotImplementedError

    def _vector_(self, R):
        raise NotImplementedError

    def _check_valid(self):
        """
        Check that this object is valid, i.e., the session in which this
        object is defined is still running. This is relevant for
        interpreters that can't be interrupted via ctrl-C, hence get
        restarted.
        """
        try:
            P = self.parent()
            if P is None:
                raise ValueError("The %s session in which this object was defined is no longer running." % P.name())
        except AttributeError:
            raise ValueError("The session in which this object was defined is no longer running.")
        return P

    def __del__(self):
        try:
            self._check_valid()
        except ValueError:
            return
        if hasattr(self, '_name'):
            P = self.parent()
            if P is not None:
                P.clear(self._name)

    def _sage_repr(self):
        """
        Return a sage-friendly string representation of the object.

        Some programs use different notation to Sage, e.g. Mathematica
        writes lists with {} instead of [].  This method calls repr(self)
        then converts the foreign notation into Sage's notation.

        OUTPUT:

        A string representation of the object that is ready for
        sage_eval().

        EXAMPLES::

            sage: repr(mathematica([1,2,3]))    # optional - mathematica
            '{1, 2, 3}'
            sage: mathematica([1,2,3])._sage_repr() # optional - mathematica
            '[1, 2, 3]'

        ::

            sage: gp(10.^80)._sage_repr()
            '1.0000000000000000000000000000000000000e80'    # 64-bit
            '1.000000000000000000000000000e80'              # 32-bit
            sage: mathematica('10.^80')._sage_repr()  # optional - mathematica
            '1.e80'

        AUTHORS:

        - Felix Lawrence (2009-08-21)
        """
        # TO DO: this could use file transfers when self.is_remote()

        string = repr(self).replace('\n', ' ').replace('\r', '')
        # Translate the external program's function notation to Sage's
        lfd = self.parent()._left_func_delim()
        if '(' != lfd:
            string = string.replace(lfd, '(')
        rfd = self.parent()._right_func_delim()
        if ')' != rfd:
            string = string.replace(rfd, ')')
        # Translate the external program's list formatting to Sage's
        lld = self.parent()._left_list_delim()
        if '[' != lld:
            string = string.replace(lld, '[')
        rld = self.parent()._right_list_delim()
        if ']' != rld:
            string = string.replace(rld, ']')
        # Translate the external program's exponent formatting
        expl = self.parent()._exponent_symbol()
        if 'e' != expl:
            string = string.replace(expl, 'e')
        return string

    def _sage_(self):
        """
        Attempt to return a Sage version of this object.
        This is a generic routine that just tries to evaluate
        the repr(self).

        EXAMPLES::

            sage: gp(1/2)._sage_()
            1/2
            sage: _.parent()
            Rational Field

        AUTHORS:

        - William Stein

        - Felix Lawrence (2009-08-21)
        """
        string = self._sage_repr()
        try:
            return sage.misc.sage_eval.sage_eval(string)
        except Exception:
            raise NotImplementedError("Unable to parse output: %s" % string)

    def sage(self, *args, **kwds):
        """
        Attempt to return a Sage version of this object.

        This method does nothing more than calling :meth:`_sage_`,
        simply forwarding any additional arguments.

        EXAMPLES::

            sage: gp(1/2).sage()
            1/2
            sage: _.parent()
            Rational Field
            sage: singular.lib("matrix")
            sage: R = singular.ring(0, '(x,y,z)', 'dp')
            sage: singular.matrix(2,2).sage()
            [0 0]
            [0 0]
        """
        return self._sage_(*args, **kwds)

    def __repr__(self):
        """
        To obtain the string representation, it is first checked whether
        the element is still valid. Then, if ``self._cached_repr`` is
        a string then it is returned. Otherwise, ``self._repr_()``
        is called (and the result is cached, if ``self._cached_repr``
        evaluates to ``True``).

        If the string obtained so far contains ``self._name``, then it
        is replaced by ``self``'s custom name, if available.

        To implement a custom string representation, override the method
        ``_repr_``, but do not override this double underscore method.

        EXAMPLES:

        Here is one example showing that the string representation will
        be cached when requested::

            sage: from sage.interfaces.maxima_lib import maxima_lib
            sage: M = maxima_lib('sqrt(2) + 1/3')
            sage: M._cached_repr
            True
            sage: repr(M) is repr(M)   # indirect doctest
            True
            sage: M._cached_repr
            'sqrt(2)+1/3'
            sage: M
            sqrt(2)+1/3

        If the interface breaks then it is reflected in the string representation::

            sage: s = singular('2')
            sage: s
            2
            sage: singular.quit()
            sage: s
            (invalid Singular object -- The singular session in which this object was defined is no longer running.)
        """
        try:
            self._check_valid()
        except ValueError as msg:
            return '(invalid {} object -- {})'.format(self.parent() or type(self), msg)
        cr = getattr(self, '_cached_repr', None)
        if isinstance(cr, str):
            s = cr
        else:
            s = self._repr_()
        if self._name in s and self.get_custom_name() is not None:
            s = s.replace(self._name, self.get_custom_name())
        if cr:
            self._cached_repr = s
        return s

    def _repr_(self):
        """
        Default implementation of a helper method for string representation.

        It is supposed that immediately before calling this method,
        the validity of ``self``'s parent was confirmed. So, when you
        override this method, you can assume that the parent is valid.

        TESTS:

        In :issue:`22501`, several string representation methods have been
        removed in favour of using the default implementation. The corresponding
        tests have been moved here::

            sage: gap(SymmetricGroup(8))    # indirect doctest
            SymmetricGroup( [ 1 .. 8 ] )
            sage: gap(2)
            2
            sage: x = var('x')
            sage: giac(x)
            sageVARx
            sage: giac(5)
            5
            sage: M = matrix(QQ,2,range(4))
            sage: giac(M)
            [[0,1],[2,3]]
            sage: x = var('x')                  # optional - maple
            sage: maple(x)                      # optional - maple
            x
            sage: maple(5)                      # optional - maple
            5
            sage: M = matrix(QQ,2,range(4))     # optional - maple
            sage: maple(M)                      # optional - maple
            Matrix(2, 2, [[0,1],[2,3]])
            sage: maxima('sqrt(2) + 1/3')
            sqrt(2)+1/3
            sage: mupad.package('"MuPAD-Combinat"')  # optional - mupad-Combinat
            sage: S = mupad.examples.SymmetricFunctions(); S # optional - mupad-Combinat
            examples::SymmetricFunctions(Dom::ExpressionField())
        """
        P = self.parent()
        try:
            if self._get_using_file:
                return P.get_using_file(self._name).rstrip()
        except AttributeError:
            return self.parent().get(self._name).rstrip()

    def __getattr__(self, attrname):
        try:
            P = self._check_valid()
        except ValueError:
            raise AttributeError(attrname)
        if attrname[:1] == "_":
            raise AttributeError
        return P._function_element_class()(self, attrname)

    def get_using_file(self):
        """
        Return this element's string representation using a file. Use this
        if ``self`` has a huge string representation. It'll be way faster.

        EXAMPLES::

            sage: a = maxima(str(2^1000))
            sage: a.get_using_file()
            '10715086071862673209484250490600018105614048117055336074437503883703510511249361224931983788156958581275946729175531468251871452856923140435984577574698574803934567774824230985421074605062371141877954182153046474983581941267398767559165543946077062914571196477686542167660429831652624386837205668069376'
        """
        try:
            self._check_valid()
        except ValueError as msg:
            return '(invalid {} object -- {})'.format(self.parent() or type(self), msg)
        return self.parent().get_using_file(self._name)

    def hasattr(self, attrname):
        """
        Return whether the given attribute is already defined by this
        object, and in particular is not dynamically generated.

        EXAMPLES::

            sage: m = maxima('2')
            sage: m.hasattr('integral')
            True
            sage: m.hasattr('gcd')
            False
        """
        return not isinstance(getattr(self, attrname), (InterfaceFunctionElement, InterfaceElement))

    def attribute(self, attrname):
        """
        If this wraps the object x in the system, this returns the object
        x.attrname. This is useful for some systems that have object
        oriented attribute access notation.

        EXAMPLES::

            sage: g = gap('SO(1,4,7)')
            sage: k = g.InvariantQuadraticForm()
            sage: k.attribute('matrix')
            [ [ 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ],
              [ 0*Z(7), 0*Z(7), Z(7), 0*Z(7) ], [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ] ]

        ::

            sage: e = gp('ellinit([0,-1,1,-10,-20])')
            sage: e.attribute('j')
            -122023936/161051
        """
        P = self._check_valid()
        return P('%s.%s' % (self.name(), attrname))

    def __getitem__(self, n):
        P = self._check_valid()
        if not isinstance(n, tuple):
            return P.new('%s[%s]' % (self._name, n))
        else:
            return P.new('%s[%s]' % (self._name, str(n)[1:-1]))

    def __int__(self):
        """
        EXAMPLES::

            sage: int(maxima('1'))
            1
            sage: type(_)
            <... 'int'>
        """
        return int(repr(self))

    def bool(self):
        """
        Convert this element to a boolean.

        EXAMPLES::

            sage: singular(0).bool()
            False
            sage: singular(1).bool()
            True
        """
        return bool(self)

    def __bool__(self):
        """
        Return whether this element is not ``False``.

        .. NOTE::

            This method needs to be overridden if the subprocess would
            not return a string representation of a boolean value unless
            an explicit print command is used.

        EXAMPLES::

            sage: bool(maxima(0))
            False
            sage: bool(maxima(1))
            True

        TESTS:

        By default this returns ``True`` for elements that are considered to be
        not ``False`` by the interface (:issue:`28705`)::

            sage: bool(giac('"a"'))
            True
        """
        P = self._check_valid()
        cmd = '%s %s %s' % (self._name, P._equality_symbol(),
                            P._false_symbol())
        return P.eval(cmd) != P._true_symbol()

    def __float__(self):
        """
        EXAMPLES::

            sage: m = maxima('1/2')
            sage: m.__float__()
            0.5
            sage: float(m)
            0.5
        """
        return float(repr(self))

    def _integer_(self, ZZ=None):
        """
        EXAMPLES::

            sage: m = maxima('1')
            sage: m._integer_()
            1
            sage: _.parent()
            Integer Ring
            sage: QQ(m)
            1
        """
        from sage.rings.integer import Integer
        return Integer(repr(self))

    def _rational_(self):
        """
        EXAMPLES::

            sage: m = maxima('1/2')
            sage: m._rational_()
            1/2
            sage: _.parent()
            Rational Field
            sage: QQ(m)
            1/2
        """
        from sage.rings.rational import Rational
        return Rational(repr(self))

    def name(self, new_name=None):
        """
        Return the name of ``self``. If ``new_name`` is passed in, then this
        function returns a new object identical to ``self`` whose name is
        ``new_name``.

        Note that this can overwrite existing variables in the system.

        EXAMPLES::

            sage: # optional - rpy2
            sage: x = r([1,2,3]); x
            [1] 1 2 3
            sage: x.name()
            'sage...'
            sage: x = r([1,2,3]).name('x'); x
            [1] 1 2 3
            sage: x.name()
            'x'

        ::

            sage: s5 = gap.SymmetricGroup(5).name('s5')
            sage: s5
            SymmetricGroup( [ 1 .. 5 ] )
            sage: s5.name()
            's5'
        """
        if new_name is not None:
            if not isinstance(new_name, str):
                raise TypeError("new_name must be a string")
            p = self.parent()
            p.set(new_name, self._name)
            return p._object_class()(p, new_name, is_name=True)

        return self._name

    def gen(self, n):
        P = self._check_valid()
        return P.new('%s.%s' % (self._name, int(n)))

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

            sage: a = gp('23')
            sage: b = gp('5')
            sage: a._operation('%', b)
            3
            sage: a._operation('19+')
            42
            sage: a._operation('!@#$%')
            Traceback (most recent call last):
            ...
            TypeError: Error executing code in GP:...
        """
        P = self._check_valid()
        if other is None:
            cmd = '%s %s' % (operation, self._name)
        else:
            cmd = '%s %s %s' % (self._name, operation, other._name)
        try:
            return P.new(cmd)
        except Exception as msg:
            raise TypeError(msg)

    def _add_(self, right):
        """
        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: g = maxima.sin(x)
            sage: f + g
            sin(_SAGE_VAR_x)+cos(_SAGE_VAR_x)
            sage: f + 2
            cos(_SAGE_VAR_x)+2
            sage: 2 + f
            cos(_SAGE_VAR_x)+2

        ::

            sage: x,y = var('x,y')
            sage: f = maxima.function('x','sin(x)')
            sage: g = maxima.function('x','-cos(x)')
            sage: f+g
            sin(x)-cos(x)
            sage: f+3
            sin(x)+3

        The Maxima variable ``x`` is different from the Sage symbolic variable::

            sage: (f+maxima.cos(x))
            sin(x)+cos(_SAGE_VAR_x)
            sage: (f+maxima.cos(y))
            sin(x)+cos(_SAGE_VAR_y)

        Note that you may get unexpected results when calling symbolic expressions
        and not explicitly giving the variables::

            sage: (f+maxima.cos(x))(2)
            cos(_SAGE_VAR_x)+sin(2)
            sage: (f+maxima.cos(y))(2)
            cos(_SAGE_VAR_y)+sin(2)
        """
        return self._operation("+", right)

    def _sub_(self, right):
        """
        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: g = maxima.sin(x)
            sage: f - g
            cos(_SAGE_VAR_x)-sin(_SAGE_VAR_x)
            sage: f - 2
            cos(_SAGE_VAR_x)-2
            sage: 2 - f
            2-cos(_SAGE_VAR_x)

        ::

            sage: x,y = var('x,y')
            sage: f = maxima.function('x','sin(x)')

        The Maxima variable ``x`` is different from the Sage symbolic variable::

            sage: (f-maxima.cos(x))
            sin(x)-cos(_SAGE_VAR_x)
            sage: (f-maxima.cos(y))
            sin(x)-cos(_SAGE_VAR_y)

        Note that you may get unexpected results when calling symbolic expressions
        and not explicitly giving the variables::

            sage: (f-maxima.cos(x))(2)
            sin(2)-cos(_SAGE_VAR_x)
            sage: (f-maxima.cos(y))(2)
            sin(2)-cos(_SAGE_VAR_y)
        """
        return self._operation('-', right)

    def _neg_(self):
        """
        EXAMPLES::

            sage: f = maxima('sin(x)')
            sage: -f
            -sin(x)
            sage: f = maxima.function('x','sin(x)')
            sage: -f
            -sin(x)
        """
        return self._operation('-')

    def _mul_(self, right):
        """
        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: g = maxima.sin(x)
            sage: f*g
            cos(_SAGE_VAR_x)*sin(_SAGE_VAR_x)
            sage: 2*f
            2*cos(_SAGE_VAR_x)

        ::

            sage: f = maxima.function('x','sin(x)')
            sage: g = maxima('cos(x)') # not a function!
            sage: f*g
            cos(x)*sin(x)
            sage: _(2)
            cos(2)*sin(2)

        ::

            sage: f = maxima.function('x','sin(x)')
            sage: g = maxima('cos(x)')
            sage: g*f
            cos(x)*sin(x)
            sage: _(2)
            cos(2)*sin(2)
            sage: 2*f
            2*sin(x)
        """
        return self._operation('*', right)

    def _div_(self, right):
        """
        EXAMPLES::

            sage: f = maxima.cos(x)
            sage: g = maxima.sin(x)
            sage: f/g
            cos(_SAGE_VAR_x)/sin(_SAGE_VAR_x)
            sage: f/2
            cos(_SAGE_VAR_x)/2

        ::

            sage: f = maxima.function('x','sin(x)')
            sage: g = maxima('cos(x)')
            sage: f/g
            sin(x)/cos(x)
            sage: _(2)
            sin(2)/cos(2)

        ::

            sage: f = maxima.function('x','sin(x)')
            sage: g = maxima('cos(x)')
            sage: g/f
            cos(x)/sin(x)
            sage: _(2)
            cos(2)/sin(2)
            sage: 2/f
            2/sin(x)
        """
        return self._operation("/", right)

    def __invert__(self):
        """
        EXAMPLES::

            sage: f = maxima('sin(x)')
            sage: ~f
            1/sin(x)
            sage: f = maxima.function('x','sin(x)')
            sage: ~f
            1/sin(x)
        """
        return self._operation('1/')

    def _mod_(self, right):
        """
        EXAMPLES::

            sage: f = gp("x^3 + x")
            sage: g = gp("2*x + 1")
            sage: f % g
            -5/8
        """
        return self._operation("%", right)

    def __pow__(self, n):
        """
        EXAMPLES::

            sage: a = maxima('2')
            sage: a^(3/4)
            2^(3/4)

        ::

            sage: f = maxima.function('x','sin(x)')
            sage: g = maxima('-cos(x)')
            sage: f^g
            1/sin(x)^cos(x)

        ::

            sage: f = maxima.function('x','sin(x)')
            sage: g = maxima('-cos(x)') # not a function
            sage: g^f
            (-cos(x))^sin(x)
        """
        P = self._check_valid()
        if parent(n) is not P:
            n = P(n)
        return self._operation("^", n)
