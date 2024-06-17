"""
Mixin For Extra Tab Completions

The :class:`ExtraTabCompletion` class helps you to extend the tab
completions for objects. This is used in interfaces with third-party
interpreters where we want, for example, ``gap.[TAB]`` to also list
GAP commands and not only the Python methods of the gap interface. To
use it, just inherit (usually by multiple inheritance) from
:class:`ExtraTabCompletion` and implement a ``_tab_completion`` method
that returns a list of strings. These strings are then included in the
tab completions on instances. It is up to you to make these "fake"
method names work, usually by implementing ``__getattr__``.

EXAMPLES::

    sage: from sage.interfaces.tab_completion import ExtraTabCompletion
    sage: class Foo(ExtraTabCompletion, object):
    ....:      def a(self):
    ....:          return 1
    ....:      def _tab_completion(self):
    ....:          return ['c', 'd']
    sage: f = Foo()
    sage: f.b = 2
    sage: sorted(dir(f))
    [..., '_tab_completion', 'a', 'b', 'c', 'd']
"""
import builtins


class ExtraTabCompletion():

    def __dir__(self):
        """
        Add to the ``dir()`` output.

        This is used by IPython to read off the tab completions.

        EXAMPLES::

            sage: from sage.interfaces.tab_completion import ExtraTabCompletion
            sage: obj = ExtraTabCompletion()
            sage: dir(obj)
            Traceback (most recent call last):
            ...
            NotImplementedError: <class 'sage.interfaces.tab_completion.ExtraTabCompletion'> must implement _tab_completion() method
        """
        try:
            tab_fn = self._tab_completion
        except AttributeError:
            raise NotImplementedError(
                '{0} must implement _tab_completion() method'.format(self.__class__))
        return dir(self.__class__) + list(self.__dict__) + tab_fn()


def completions(s, globs):
    """
    Return a list of completions in the given context.

    INPUT:

    - ``s`` -- string

    - ``globs`` -- string: object dictionary; context in which to
      search for completions, e.g., :func:`globals()`

    OUTPUT: list of strings

    EXAMPLES::

         sage: X.<x> = PolynomialRing(QQ)
         sage: import sage.interfaces.tab_completion as s
         sage: p = x**2 + 1
         sage: s.completions('p.co',globals()) # indirect doctest
         ['p.coefficients',...]

         sage: s.completions('dic',globals()) # indirect doctest
         ['dickman_rho', 'dict']
    """
    if not s:
        raise ValueError('empty string')

    if '.' not in s:
        n = len(s)
        v = [x for x in globs if x[:n] == s]
        v += [x for x in builtins.__dict__ if x[:n] == s]
    else:
        i = s.rfind('.')
        method = s[i + 1:]
        obj = s[:i]
        n = len(method)
        try:
            O = eval(obj, globs)
            D = dir(O)
            if not method:
                v = [obj + '.' + x for x in D if x and x[0] != '_']
            else:
                v = [obj + '.' + x for x in D if x[:n] == method]
        except Exception:
            v = []
    return sorted(set(v))
