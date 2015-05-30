# -*- coding: utf-8 -*-
"""
Lazy strings

Based on speaklater: https://github.com/mitsuhiko/speaklater.

A lazy string is an object that behaves almost exactly like a string
but where the value is not computed until needed.  To define a lazy
string you specify a function that produces a string together with the
appropriate arguments for that function.  Sage uses lazy strings in
:mod:`sage.misc.misc` so that the filenames for SAGE_TMP (which
depends on the pid of the process running Sage) are not computed when
importing the Sage library.  This means that when the doctesting code
imports the Sage library and then forks, the variable SAGE_TMP depends
on the new pid rather than the old one.

EXAMPLES::

    sage: from sage.misc.lazy_string import lazy_string
    sage: L = []
    sage: s = lazy_string(lambda x: str(len(x)), L)
    sage: L.append(5)
    sage: s
    l'1'

Note that the function is recomputed each time::

    sage: L.append(6)
    sage: s
    l'2'
"""

#Copyright (c) 2009 by Armin Ronacher.
#
#Some rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are
#met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above
#      copyright notice, this list of conditions and the following
#      disclaimer in the documentation and/or other materials provided
#      with the distribution.
#
#    * The names of the contributors may not be used to endorse or
#      promote products derived from this software without specific
#      prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import types

def is_lazy_string(obj):
    """
    Checks if the given object is a lazy string.

    EXAMPLES::

        sage: from sage.misc.lazy_string import lazy_string, is_lazy_string
        sage: f = lambda: "laziness"
        sage: s = lazy_string(f)
        sage: is_lazy_string(s)
        True
    """
    return isinstance(obj, _LazyString)

def lazy_string(func, *args, **kwargs):
    """
    Creates a lazy string by invoking func with args.

    EXAMPLES::

        sage: from sage.misc.lazy_string import lazy_string
        sage: f = lambda: "laziness"
        sage: s = lazy_string(f); s
        l'laziness'
    """
    return _LazyString(func, args, kwargs)

def _make_lazy_string(ftype, fpickle, args, kwargs):
    """
    Used for pickling.

    EXAMPLES::

        sage: from sage.misc.lazy_string import _make_lazy_string
        sage: s = _make_lazy_string(None, lambda: "laziness", (), {})
        sage: s
        l'laziness'
    """
    if ftype == 'func':
        from sage.misc.fpickle import unpickle_function
        f = unpickle_function(fpickle)
    else:
        f = fpickle
    return _LazyString(f, args, kwargs)

class _LazyString(object):
    """
    Class for strings created by a function call.

    The proxy implementation attempts to be as complete as possible, so that
    the lazy objects should mostly work as expected, for example for sorting.

    EXAMPLES::

        sage: from sage.misc.lazy_string import lazy_string
        sage: f = lambda x: "laziness in the " + repr(x)
        sage: s = lazy_string(f, ZZ); s
        l'laziness in the Integer Ring'
    """
    __slots__ = ('_func', '_args', '_kwargs')

    def __init__(self, func, args, kwargs):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda x: "laziness" + repr(x)
            sage: s = lazy_string(f, 5); s
            l'laziness5'
        """
        self._func = func
        self._args = args
        self._kwargs = kwargs

    value = property(lambda x: x._func(*x._args, **x._kwargs))

    def __contains__(self, key):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: 'zi' in s
            True
            sage: 'ni' in s
            False
        """
        return key in self.value

    def __nonzero__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: bool(lazy_string(f))
            True
            sage: f = lambda: ""
            sage: bool(lazy_string(f))
            False
        """
        return bool(self.value)

    def __dir__(self):
        """
        We assume that the underlying value provides the methods of a
        unicode string.

        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: "split" in dir(s) # indirect doctest
            True
        """
        return dir(unicode)

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: "".join(list(s)) # indirect doctest
            'laziness'
        """
        return iter(self.value)

    def __len__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: len(s)
            8
        """
        return len(self.value)

    def __str__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: str(s) # indirect doctest
            'laziness'
        """
        return str(self.value)

    def __unicode__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: unicode(s) # indirect doctest
            u'laziness'
        """
        return unicode(self.value)

    def __add__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s + " supreme"
            'laziness supreme'
        """
        return self.value + other

    def __radd__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: "no " + s
            'no laziness'
        """
        return other + self.value

    def __mod__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laz%sss"
            sage: s = lazy_string(f)
            sage: s % "ine"
            'laziness'
        """
        return self.value % other

    def __rmod__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "ine"
            sage: s = lazy_string(f)
            sage: "laz%sss" % s
            'laziness'
        """
        return other % self.value

    def __mul__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s * 2
            'lazinesslaziness'
        """
        return self.value * other

    def __rmul__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: 2 * s
            'lazinesslaziness'
        """
        return other * self.value

    def __lt__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s < 'laziness'
            False
            sage: s < 'azi'
            False
            sage: s < s
            False
        """
        return self.value < other

    def __le__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s <= 'laziness'
            True
            sage: s <= 'azi'
            False
            sage: s <= s
            True
        """
        return self.value <= other

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s == 'laziness'
            True
            sage: s == 'azi'
            False
            sage: s == s
            True
        """
        return self.value == other

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s != 'laziness'
            False
            sage: s != 'azi'
            True
            sage: s != s
            False
        """
        return self.value != other

    def __gt__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s > 'laziness'
            False
            sage: s > 'azi'
            True
            sage: s > s
            False
        """
        return self.value > other

    def __ge__(self, other):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s >= 'laziness'
            True
            sage: s >= 'azi'
            True
            sage: s >= s
            True
        """
        return self.value >= other

    def __getattr__(self, name):
        """
        We pass attribute lookup through to the underlying value.

        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s.endswith('ess')
            True
            sage: s.find('i')
            3
        """
        if name == '__members__':
            return self.__dir__()
        return getattr(self.value, name)

    def __reduce__(self):
        """
        Pickling.

        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: TestSuite(s).run() # indirect doctest
        """
        import types
        if isinstance(self._func, types.FunctionType):
            from sage.misc.fpickle import pickle_function
            f = pickle_function(self._func)
            ftype = 'func'
        else:
            f = self.func
            ftype = None
        return _make_lazy_string, (ftype, f, self._args, self._kwargs)

    def __getitem__(self, key):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s[4]
            'n'
        """
        return self.value[key]

    def __copy__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: copy(s) is s
            True
        """
        return self

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.misc.lazy_string import lazy_string
            sage: f = lambda: "laziness"
            sage: s = lazy_string(f)
            sage: s # indirect doctest
            l'laziness'
        """
        try:
            return 'l' + repr(self.value)
        except Exception:
            return '<%s broken>' % self.__class__.__name__
