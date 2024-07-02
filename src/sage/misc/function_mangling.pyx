# sage_setup: distribution = sagemath-objects
# Copyright (c) 2009, Tom Boothby <boothby@math.washington.edu>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of Sage nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY Tom Boothby ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL Tom Boothby BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

Function Mangling

This module provides utilities for extracting information about python
functions.

AUTHORS:

- Tom Boothby (2009): Original version in Python
- Simon King (2011): Use Cython. Speedup of ``fix_to_pos``, cleaning documentation.
"""

from sage.misc.sageinspect import sage_getargspec

cdef class ArgumentFixer:
    """
    This class provides functionality to normalize the arguments
    passed into a function.  While the various ways of calling a
    function are perfectly equivalent from the perspective of the
    callee, they don't always look the same for an object
    watching the caller.  For example,
    ::

        sage: def f(x = 10):
        ....:     return min(1,x)

    the following calls are equivalent,
    ::

        sage: f()
        1
        sage: f(10)
        1
        sage: f(x=10)
        1

    but from the perspective of a wrapper, they are different::

        sage: def wrap(g):
        ....:     def _g(*args,**kwargs):
        ....:         print("{} {}".format(args, kwargs))
        ....:         return g(*args, **kwargs)
        ....:     return _g
        sage: h = wrap(f)
        sage: t = h()
        () {}
        sage: t = h(10)
        (10,) {}
        sage: t = h(x=10)
        () {'x': 10}

    For the purpose of cached functions, it is important not
    to distinguish between these uses.

    INPUT:

    - ``f`` -- a function
    - ``classmethod`` -- boolean (default: ``False``); ``True`` if the function
      is a classmethod and therefore the first argument is expected to be the
      class instance. In that case, we ignore the first argument.

    EXAMPLES::

        sage: from sage.misc.function_mangling import ArgumentFixer
        sage: def wrap2(g):
        ....:     af = ArgumentFixer(g)
        ....:     def _g(*args, **kwargs):
        ....:         print(af.fix_to_pos())
        ....:         return g(*args,**kwargs)
        ....:     return _g
        sage: h2 = wrap2(f)
        sage: t = h2()
        ((10,), ())
        sage: t = h2(10)
        ((10,), ())
        sage: t = h2(x=10)
        ((10,), ())

    ::

        sage: class one:
        ....:    def __init__(self, x = 1):
        ....:       self.x = x
        sage: af = ArgumentFixer(one.__init__, classmethod=True)
        sage: af.fix_to_pos(1,2,3,a=31,b=2,n=3)
        ((1, 2, 3), (('a', 31), ('b', 2), ('n', 3)))
    """
    def __init__(self, f, classmethod = False):
        try:
            arg_names, varargs, varkw, defaults, kwonlyargs, kwonlydefaults, annotations = sage_getargspec(f)
        except AttributeError:
            # This error occurs if f is defined in a Cython file and the
            # source file has gone.
            if classmethod:
                arg_names = ['self']
                varargs = 'args'
                varkws = 'kwds'
                defaults = None
            else:
                arg_names = []
                varargs = 'args'
                varkws = 'kwds'
                defaults = None
        if defaults is None:
            self._default_tuple = defaults = ()
        else:
            self._default_tuple = tuple(defaults)

        #code = f.__code__

        self.f = f
        self._ndefault = len(defaults)
        if classmethod:
            self._nargs = len(arg_names)-1 #code.co_argcount-1
            self._arg_names = tuple(arg_names[1:]) #code.co_varnames[1:self._nargs+1]
        else:
            self._nargs = len(arg_names) #code.co_argcount
            self._arg_names = tuple(arg_names) #code.co_varnames[:self._nargs]
        self._classmethod = classmethod

        cdef dict default_map
        self._defaults = default_map = {}
        for k,v in zip(self._arg_names[-self._ndefault:], defaults):
            default_map[k] = v

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.misc.function_mangling import ArgumentFixer
            sage: g = ArgumentFixer(number_of_partitions); g                            # needs sage.combinat
            Argument Fixer of <function number_of_partitions at 0x...>
        """
        return "Argument Fixer of %s"%self.f

    def fix_to_named(self, *args,**kwargs):
        """
        Normalize the arguments with a preference for named arguments.

        INPUT:

        - any positional and named arguments.

        OUTPUT: we return a tuple

            `(e_1, e_2, ..., e_k), ((n_1, v_1), ... , (n_m, v_m))`

        where `n_1, ... , n_m` are the names of the arguments and
        `v_1, ..., v_m` are the values passed in; and `e_1, ..., e_k` are
        the unnamed arguments.  We minimize `k`.

        The defaults are extracted from the function and filled
        into the list ``K`` of named arguments. The names `n_1, ..., n_t`
        are in order of the function definition, where `t` is the number
        of named arguments.  The remaining names, `n_{t+1}, ..., n_m` are
        given in alphabetical order.  This is useful to extract
        the names of arguments, but **does not** maintain
        equivalence of
        ::

            A,K = self.fix_to_pos(...)
            self.f(*A,**dict(K))`

        and
        ::

            self.f(...)

        in all cases.

        EXAMPLES::

            sage: from sage.misc.function_mangling import ArgumentFixer
            sage: def sum3(a, b, c=3, *args, **kwargs):
            ....:     return a + b + c
            sage: AF = ArgumentFixer(sum3)
            sage: AF.fix_to_named(1, 2, 3, 4, 5, 6, f=14, e=16)
            ((4, 5, 6), (('a', 1), ('b', 2), ('c', 3), ('e', 16), ('f', 14)))
            sage: AF.fix_to_named(1,2,f=14)
            ((), (('a', 1), ('b', 2), ('c', 3), ('f', 14)))
        """
        cdef list ARGS = []
        cdef tuple arg_names = self._arg_names
        cdef int lenargs = len(args)
        cdef dict defaults = self._defaults
        cdef int i
        cdef dict kwargs_ = dict(kwargs)
        for i from 0<=i<self._nargs:
            name = arg_names[i]
            if i >= lenargs:
                if name in kwargs_:
                    val = kwargs_[name]
                    del kwargs_[name]
                else:
                    val = defaults[name]
            else:
                val = args[i]
            ARGS.append((name,val))
        extra_args = args[self._nargs:]
        for k in sorted(kwargs_.keys()):
            ARGS.append((k,kwargs_[k]))
        return tuple(extra_args), tuple(ARGS)

    def fix_to_pos(self, *args, **kwds):
        """
        Normalize the arguments with a preference for positional arguments.

        INPUT:

        - ``*args``, ``**kwds`` -- any positional or named arguments

        OUTPUT: we return a tuple

            `(e_1, e_2, ..., e_k), ((n_1, v_1), ... , (n_m, v_m))`

        where `n_1, ... , n_m` are the names of the arguments and
        `v_1, ..., v_m` are the values passed in; and `e_1, ..., e_k`
        are the unnamed arguments. We minimize `m`.

        The commands
        ::

            A,K = self.fix_to_pos(...)
            self.f(*A,**dict(K))

        are equivalent to
        ::

            self.f(...)

        though defaults are extracted from the function and
        appended to the tuple ``A`` of positional arguments.
        The names `n_1, ..., n_m` are given in alphabetical
        order.

        EXAMPLES::

            sage: from sage.misc.function_mangling import ArgumentFixer
            sage: def do_something(a, b, c=3, *args, **kwargs):
            ....:     print("{} {} {} {} {}".format(a, b, c, args,
            ....:                                   sorted(kwargs.items())))
            sage: AF = ArgumentFixer(do_something)
            sage: A, K = AF.fix_to_pos(1, 2, 3, 4, 5, 6, f=14, e=16)
            sage: print("{} {}".format(A, K))
            (1, 2, 3, 4, 5, 6) (('e', 16), ('f', 14))
            sage: do_something(*A, **dict(K))
            1 2 3 (4, 5, 6) [('e', 16), ('f', 14)]
            sage: do_something(1, 2, 3, 4, 5, 6, f=14, e=16)
            1 2 3 (4, 5, 6) [('e', 16), ('f', 14)]
        """
        return self.fix_to_pos_args_kwds(args, kwds)

    cdef fix_to_pos_args_kwds(self, tuple args, dict kwds):
        """
        Fast Cython implementation of :meth:`fix_to_pos`.
        """
        cdef Py_ssize_t lenargs = len(args)
        cdef Py_ssize_t nargs = self._nargs
        cdef tuple arg_names = self._arg_names
        cdef dict defaults = self._defaults
        # a shortpath for the case of no named arguments:
        if not kwds:
            if lenargs >= nargs:
                return args, ()
            # we take the given arguments, plus the default arguments
            return args + self._default_tuple[-nargs+lenargs:],()
        cdef list Largs = list(args)
        cdef dict kwargs = dict(kwds)
        cdef Py_ssize_t i
        for i in range(lenargs, nargs):
            # in addition to the positional arguments, we take the
            # ones with default values, unless they are overridded by
            # the named arguments.
            name = arg_names[i]
            if name in kwargs:
                val = kwargs[name]
                del kwargs[name]
            else:
                val = defaults[name]
            Largs.append(val)
        cdef list Items = sorted(kwargs.items())
        return tuple(Largs), tuple(Items)
