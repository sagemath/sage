# -*- encoding: utf-8 -*-
# cython: old_style_globals=True
r"""
Abstract base class for Sage objects
"""

from __future__ import absolute_import, print_function

from six.moves import cPickle
import os
import sys
from six.moves import cStringIO as StringIO
from sage.misc.sage_unittest import TestSuite

sys_modules = sys.modules

# change to import zlib to use zlib instead; but this
# slows down loading any data stored in the other format
import zlib; comp = zlib
import bz2; comp_other = bz2

op_LT = Py_LT   # operator <
op_LE = Py_LE   # operator <=
op_EQ = Py_EQ   # operator ==
op_NE = Py_NE   # operator !=
op_GT = Py_GT   # operator >
op_GE = Py_GE   # operator >=


from sage.misc.superseded import deprecated_function_alias
py_rich_to_bool = deprecated_function_alias(21128, rich_to_bool)


cdef process(s):
    if s[-5:] != '.sobj':
        return s + '.sobj'
    else:
        return s


cdef class SageObject:
    """
    Base class for all (user-visible) objects in Sage

    Every object that can end up being returned to the user should
    inherit from :class:`SageObject`.

    .. automethod:: _ascii_art_
    .. automethod:: _cache_key
    """
    def _test_new(self, **options):
        """
        Check that ``cls.__new__(cls)`` does not crash Python,
        where ``cls = type(self)``.

        It is perfectly legal for ``__new__`` to raise ordinary
        exceptions.

        EXAMPLES::

            sage: SageObject()._test_new()
        """
        cdef type cls = type(self)
        try:
            cls.__new__(cls)
        except Exception:
            pass

    #######################################################################
    # Textual representation code
    #######################################################################

    def rename(self, x=None):
        r"""
        Change self so it prints as x, where x is a string.

        .. NOTE::

           This is *only* supported for Python classes that derive
           from SageObject.

        EXAMPLES::

            sage: x = PolynomialRing(QQ, 'x', sparse=True).gen()
            sage: g = x^3 + x - 5
            sage: g
            x^3 + x - 5
            sage: g.rename('a polynomial')
            sage: g
            a polynomial
            sage: g + x
            x^3 + 2*x - 5
            sage: h = g^100
            sage: str(h)[:20]
            'x^300 + 100*x^298 - '
            sage: h.rename('x^300 + ...')
            sage: h
            x^300 + ...

        Real numbers are not Python classes, so rename is not supported::

            sage: a = 3.14
            sage: type(a)
            <type 'sage.rings.real_mpfr.RealLiteral'>
            sage: a.rename('pi')
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support renaming: 3.14000000000000

        .. NOTE::

           The reason C-extension types are not supported by default
           is if they were then every single one would have to carry
           around an extra attribute, which would be slower and waste
           a lot of memory.

           To support them for a specific class, add a
           ``cdef public __custom_name`` attribute.
        """
        if x is None:
            #if hasattr(self, '__custom_name'):
            # that's tested in reset_name anyway...
            self.reset_name()
        else:
            try:
                self.__custom_name = str(x)
            except AttributeError:
                raise NotImplementedError("object does not support renaming: %s" % self)

    def reset_name(self):
        """
        Remove the custom name of an object.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: P
            Univariate Polynomial Ring in x over Rational Field
            sage: P.rename('A polynomial ring')
            sage: P
            A polynomial ring
            sage: P.reset_name()
            sage: P
            Univariate Polynomial Ring in x over Rational Field
        """
        if hasattr(self, '__custom_name'):
            del self.__custom_name


    def __repr__(self):
        """
        Default method for string representation.

        .. NOTE::

            Do not overwrite this method. Instead, implement
            a ``_repr_`` (single underscore) method.

        EXAMPLES:

        By default, the string representation coincides with
        the output of the single underscore ``_repr_``::

            sage: P.<x> = QQ[]
            sage: repr(P) == P._repr_()  #indirect doctest
            True

        Using :meth:`rename`, the string representation can
        be customized::

            sage: P.rename('A polynomial ring')
            sage: repr(P) == P._repr_()
            False

        The original behaviour is restored with :meth:`reset_name`.::

            sage: P.reset_name()
            sage: repr(P) == P._repr_()
            True
        """
        try:
            name = self.__custom_name
            if name is not None:
                return name
        except AttributeError:
            pass
        try:
            repr_func = self._repr_
        except AttributeError:
            return str(type(self))
        else:
            result = repr_func()
            if isinstance(result, unicode):
                # Py3 compatibility: allow _repr_ to return unicode
                return result.encode('utf-8')
            else:
                return result

    def _ascii_art_(self):
        r"""
        Return an ASCII art representation.

        To implement multi-line ASCII art output in a derived class
        you must override this method. Unlike :meth:`_repr_`, which is
        sometimes used for the hash key, the output of
        :meth:`_ascii_art_` may depend on settings and is allowed to
        change during runtime.

        OUTPUT:

        An :class:`~sage.typeset.ascii_art.AsciiArt` object, see
        :mod:`sage.typeset.ascii_art` for details.

        EXAMPLES:

        You can use the :func:`~sage.typeset.ascii_art.ascii_art` function
        to get the ASCII art representation of any object in Sage::

            sage: ascii_art(integral(exp(x+x^2)/(x+1), x))
              /
             |
             |   2
             |  x  + x
             | e
             | ------- dx
             |  x + 1
             |
            /

        Alternatively, you can use the ``%display ascii_art/simple`` magic to
        switch all output to ASCII art and back::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('tab = StandardTableaux(3)[2]; tab')
            [[1, 2], [3]]
            sage: shell.run_cell('%display ascii_art')
            sage: shell.run_cell('tab')
            1  2
            3
            sage: shell.run_cell('Tableaux.options(ascii_art="table", convention="French")')
            sage: shell.run_cell('tab')
            +---+
            | 3 |
            +---+---+
            | 1 | 2 |
            +---+---+
            sage: shell.run_cell('%display plain')
            sage: shell.run_cell('Tableaux.options._reset()')
            sage: shell.quit()

        TESTS::

            sage: 1._ascii_art_()
            1
            sage: type(_)
            <class 'sage.typeset.ascii_art.AsciiArt'>
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(repr(self).splitlines())

    def _unicode_art_(self):
        r"""
        Return a unicode art representation.

        To implement multi-line unicode art output in a derived class
        you must override this method. Unlike :meth:`_repr_`, which is
        sometimes used for the hash key, the output of
        :meth:`_unicode_art_` may depend on settings and is allowed to
        change during runtime.

        OUTPUT:

        An :class:`~sage.typeset.unicode_art.UnicodeArt` object, see
        :mod:`sage.typeset.unicode_art` for details.

        EXAMPLES:

        You can use the :func:`~sage.typeset.unicode_art.unicode_art` function
        to get the ASCII art representation of any object in Sage::

            sage: unicode_art(integral(exp(x+x^2)/(x+1), x))
            ⌠
            ⎮   2
            ⎮  x  + x
            ⎮ ℯ
            ⎮ ─────── dx
            ⎮  x + 1
            ⌡


        Alternatively, you can use the ``%display ascii_art/simple`` magic to
        switch all output to ASCII art and back::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('tab = StandardTableaux(3)[2]; tab')
            [[1, 2], [3]]
            sage: shell.run_cell('%display ascii_art')
            sage: shell.run_cell('tab')
            1  2
            3
            sage: shell.run_cell('Tableaux.options(ascii_art="table", convention="French")')
            sage: shell.run_cell('tab')
            +---+
            | 3 |
            +---+---+
            | 1 | 2 |
            +---+---+
            sage: shell.run_cell('%display plain')
            sage: shell.run_cell('Tableaux.options._reset()')
            sage: shell.quit()

        TESTS::

            sage: 1._unicode_art_()
            1
            sage: type(_)
            <class 'sage.typeset.unicode_art.UnicodeArt'>
        """
        ascii_art = self._ascii_art_()
        lines = map(unicode, ascii_art)
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(lines)

    def __hash__(self):
        r"""
        Not implemented: mutable objects inherit from this class

        EXAMPLES::

            sage: hash(SageObject())
            Traceback (most recent call last):
            ...
            TypeError: <type 'sage.structure.sage_object.SageObject'> is not hashable
        """
        raise TypeError("{} is not hashable".format(type(self)))

    def _cache_key(self):
        r"""
        Return a hashable key which identifies this objects for caching. The
        output must be hashable itself, or a tuple of objects which are
        hashable or define a ``_cache_key``.

        This method will only be called if the object itself is not hashable.

        Some immutable objects (such as `p`-adic numbers) cannot implement a
        reasonable hash function because their ``==`` operator has been
        modified to return ``True`` for objects which might behave differently
        in some computations::

            sage: K.<a> = Qq(9)
            sage: b = a + O(3)
            sage: c = a + 3
            sage: b
            a + O(3)
            sage: c
            a + 3 + O(3^20)
            sage: b == c
            True
            sage: b == a
            True
            sage: c == a
            False

        If such objects defined a non-trivial hash function, this would break
        caching in many places. However, such objects should still be usable in
        caches. This can be achieved by defining an appropriate
        ``_cache_key``::

            sage: hash(b)
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'sage.rings.padics.qadic_flint_CR.qAdicCappedRelativeElement'
            sage: @cached_method
            ....: def f(x): return x==a
            sage: f(b)
            True
            sage: f(c) # if b and c were hashable, this would return True
            False

            sage: b._cache_key()
            (..., ((0, 1),), 0, 1)
            sage: c._cache_key()
            (..., ((0, 1), (1,)), 0, 20)

        An implementation must make sure that for elements ``a`` and ``b``,
        if ``a != b``, then also ``a._cache_key() != b._cache_key()``.
        In practice this means that the ``_cache_key`` should always include
        the parent as its first argument::

            sage: S.<a> = Qq(4)
            sage: d = a + O(2)
            sage: b._cache_key() == d._cache_key() # this would be True if the parents were not included
            False

        """
        try:
            hash(self)
        except TypeError:
            raise TypeError("{} is not hashable and does not implement _cache_key()".format(type(self)))
        else:
            assert False, "_cache_key() must not be called for hashable elements"

    ##########################################################################
    # DATABASE Related code
    ##########################################################################

    def save(self, filename=None, compress=True):
        """
        Save self to the given filename.

        EXAMPLES::

            sage: f = x^3 + 5
            sage: f.save(os.path.join(SAGE_TMP, 'file'))
            sage: load(os.path.join(SAGE_TMP, 'file.sobj'))
            x^3 + 5
        """
        if filename is None:
            try:
                filename = self._default_filename
            except AttributeError:
                raise RuntimeError("no default filename, so it must be specified")
        filename = process(filename)
        try:
            self._default_filename = filename
        except AttributeError:
            pass
        open(filename, 'wb').write(self.dumps(compress))

    def dump(self, filename, compress=True):
        """
        Same as self.save(filename, compress)
        """
        return self.save(filename, compress=compress)

    def dumps(self, compress=True):
        r"""
        Dump ``self`` to a string ``s``, which can later be reconstituted
        as ``self`` using ``loads(s)``.

        There is an optional boolean argument ``compress`` which defaults to ``True``.

        EXAMPLES::

            sage: O=SageObject(); O.dumps()
            'x\x9ck`J.NLO\xd5+.)*M.)-\x02\xb2\x80\xdc\xf8\xfc\xa4\xac\xd4\xe4\x12\xae` \xdb\x1f\xc2,d\xd4l,d\xd2\x03\x00\xb7X\x10\xf1'
            sage: O.dumps(compress=False)
            '\x80\x02csage.structure.sage_object\nSageObject\nq\x01)\x81q\x02.'
        """
        # the protocol=2 is very important -- this enables
        # saving extensions classes (with no attributes).
        s = cPickle.dumps(self, protocol=2)
        if compress:
            return comp.compress(s)
        else:
            return s

    def db(self, name, compress=True):
        r"""
        Dumps self into the Sage database.  Use db(name) by itself to
        reload.

        The database directory is ``$HOME/.sage/db``

        TESTS::

            sage: SageObject().db("Test")
            doctest:... DeprecationWarning: db() is deprecated.
            See http://trac.sagemath.org/2536 for details.
        """
        from sage.misc.superseded import deprecation
        deprecation(2536, 'db() is deprecated.')
        from sage.misc.all import SAGE_DB
        return self.dump('%s/%s'%(SAGE_DB,name), compress=compress)

    #############################################################################
    # Category theory / structure
    #############################################################################

    def category(self):
        from sage.categories.all import Objects
        return Objects()

    def _test_category(self, **options):
        """
        Run generic tests on the method :meth:`.category`.

        See also: :class:`TestSuite`.

        EXAMPLES::

            sage: O = SageObject()
            sage: O._test_category()

        Let us now write a broken :meth:`.category` method::

            sage: class CCls(SageObject):
            ....:     def category(self):
            ....:         return 3
            sage: CC = CCls()
            sage: CC._test_category()
            Traceback (most recent call last):
            ...
            AssertionError: False is not true
        """
        from sage.categories.category import Category
        from sage.categories.objects import Objects
        tester = self._tester(**options)
        category = self.category()
        tester.assert_(isinstance(category, Category))
        tester.assert_(category.is_subcategory(Objects()))
        tester.assert_(self in category)

    def parent(self):
        """
        Return the type of ``self`` to support the coercion framework.

        EXAMPLES::

            sage: t = log(sqrt(2) - 1) + log(sqrt(2) + 1); t
            log(sqrt(2) + 1) + log(sqrt(2) - 1)
            sage: u = t.maxima_methods()
            sage: u.parent()
            <class 'sage.symbolic.maxima_wrapper.MaximaWrapper'>
        """
        return type(self)


    #############################################################################
    # Test framework
    #############################################################################

    def _tester(self, **options):
        """
        Returns a gadget attached to ``self`` providing testing utilities.

        This is used by :class:`sage.misc.sage_unittest.TestSuite` and the
        ``_test_*`` methods.

        EXAMPLES::

            sage: tester = ZZ._tester()

            sage: tester.assert_(1 == 1)
            sage: tester.assert_(1 == 0)
            Traceback (most recent call last):
            ...
            AssertionError: False is not true
            sage: tester.assert_(1 == 0, "this is expected to fail")
            Traceback (most recent call last):
            ...
            AssertionError: this is expected to fail

            sage: tester.assertEquals(1, 1)
            sage: tester.assertEquals(1, 0)
            Traceback (most recent call last):
            ...
            AssertionError: 1 != 0

        The available assertion testing facilities are the same as in
        :class:`unittest.TestCase`, which see (actually, by a slight
        abuse, tester is currently an instance of this class).

        TESTS::

            sage: ZZ._tester(tester = tester) is tester
            True
        """
        from sage.misc.sage_unittest import instance_tester
        return instance_tester(self, **options)

    def _test_not_implemented_methods(self, **options):
        """
        Checks that all required methods for this object are implemented

        TESTS::

            sage: class Abstract(SageObject):
            ....:     @abstract_method
            ....:     def bla(self):
            ....:         "returns bla"
            sage: class Concrete(Abstract):
            ....:     def bla(self):
            ....:         return 1
            sage: class IncompleteConcrete(Abstract):
            ....:     pass
            sage: Concrete()._test_not_implemented_methods()
            sage: IncompleteConcrete()._test_not_implemented_methods()
            Traceback (most recent call last):
            ...
            AssertionError: Not implemented method: bla

        """
        tester = self._tester(**options)
        try:
            # Disable warnings for the duration of the test
            import warnings
            warnings.filterwarnings('ignore')
            for name in dir(self):
                try:
                    getattr(self, name)
                except NotImplementedError:
                    # It would be best to make sure that this NotImplementedError was triggered by AbstractMethod
                    tester.fail("Not implemented method: %s"%name)
                except Exception:
                    pass
        finally: 
            # Restore warnings
            warnings.filters.pop(0)

    def _test_pickling(self, **options):
        """
        Checks that this object can be pickled and unpickled properly.

        EXAMPLES::

            sage: ZZ._test_pickling()

        .. SEEALSO::

            :func:`dumps`, :func:`loads`

        TESTS::

            sage: class Bla(SageObject): pass
            sage: Bla()._test_pickling()
            Traceback (most recent call last):
            ...
            PicklingError: Can't pickle <class '__main__.Bla'>: attribute lookup __main__.Bla failed

        TODO: for a stronger test, this could send the object to a
        remote Sage session, and get it back.
        """
        tester = self._tester(**options)
        from sage.misc.all import loads, dumps
        tester.assertEqual(loads(dumps(self)), self)

    #############################################################################
    # Coercions to interface objects
    #############################################################################

    # Sage
    def _sage_(self):
        return self

    def _interface_(self, I):
        """
        Return coercion of self to an object of the interface I.

        The result of coercion is cached, unless self is not a C
        extension class or ``self._interface_is_cached_()`` returns
        False.
        """
        c = self._interface_is_cached_()
        if c:
            try:
                X = self.__interface[I]
                X._check_valid()
                return X
            except (AttributeError, TypeError):
                try:
                    self.__interface = {}
                except AttributeError:
                    # do this because C-extension classes won't have
                    # an __interface attribute.
                    pass
            except (KeyError, ValueError):
                pass
        nm = I.name()
        init_func = getattr(self, '_%s_init_' % nm, None)
        if init_func is not None:
            s = init_func()
        else:
            try:
                s = self._interface_init_(I)
            except Exception:
                raise NotImplementedError("coercion of object %s to %s not implemented:\n%s\n%s" % (repr(self), I))
        X = I(s)
        if c:
            try:
                self.__interface[I] = X
            except AttributeError:
                pass
        return X

    def _interface_init_(self, I=None):
        return repr(self)

    def _interface_is_cached_(self):
        """
        Return True if the interface objects are cached.

        If you have an object x and do gp(x), the result is cached if
        this function returns True.
        """
        return True

    def _gap_(self, G=None):
        if G is None:
            import sage.interfaces.gap
            G = sage.interfaces.gap.gap
        return self._interface_(G)

    def _gap_init_(self):
        import sage.interfaces.gap
        I = sage.interfaces.gap.gap
        return self._interface_init_(I)

    def _gp_(self, G=None):
        if G is None:
            import sage.interfaces.gp
            G = sage.interfaces.gp.gp
        return self._interface_(G)

    def _gp_init_(self):
        return self._pari_init_()

    def _kash_(self, G=None):
        if G is None:
            import sage.interfaces.kash
            G = sage.interfaces.kash.kash
        return self._interface_(G)

    def _kash_init_(self):
        import sage.interfaces.kash
        I = sage.interfaces.kash.kash
        return self._interface_init_(I)

    def _axiom_(self, G=None):
        if G is None:
            import sage.interfaces.axiom
            G = sage.interfaces.axiom.axiom
        return self._interface_(G)

    def _axiom_init_(self):
        import sage.interfaces.axiom
        I = sage.interfaces.axiom.axiom
        return self._interface_init_(I)

    def _fricas_(self, G=None):
        if G is None:
            import sage.interfaces.fricas
            G = sage.interfaces.fricas.fricas
        return self._interface_(G)

    def _fricas_init_(self):
        import sage.interfaces.fricas
        I = sage.interfaces.fricas.fricas
        return self._interface_init_(I)

    def _giac_(self, G=None):
        if G is None:
            import sage.interfaces.giac
            G = sage.interfaces.giac.giac
        return self._interface_(G)

    def _giac_init_(self):
        import sage.interfaces.giac
        I = sage.interfaces.giac.giac
        return self._interface_init_(I)

    def _maxima_(self, G=None):
        if G is None:
            import sage.interfaces.maxima
            G = sage.interfaces.maxima.maxima
        return self._interface_(G)

    def _maxima_init_(self):
        import sage.interfaces.maxima
        I = sage.interfaces.maxima.maxima
        return self._interface_init_(I)

    def _maxima_lib_(self, G=None):
        from sage.interfaces.maxima_lib import maxima_lib
        return self._interface_(maxima_lib)

    def _maxima_lib_init_(self):
        return self._maxima_init_()

    def _magma_init_(self, magma):
        """
        Given a Magma interpreter M, return a string that evaluates in
        that interpreter to the Magma object corresponding to self.
        This function may call the magma interpreter when it runs.

        INPUT:

        - ``magma`` -- a Magma interface

        OUTPUT:

        - string

        EXAMPLES::

            sage: n = -3/7
            sage: n._magma_init_(magma)
            '-3/7'

        Some other examples that illustrate conversion to Magma.
        ::

            sage: n = -3/7
            sage: m2 = Magma()
            sage: magma(n)                        # optional - magma
            -3/7
            sage: magma(n).parent()               # optional - magma
            Magma
            sage: magma(n).parent() is m2         # optional - magma
            False
            sage: magma(n).parent() is magma      # optional - magma
            True

        This example illustrates caching, which happens automatically
        since K is a Python object::

            sage: K.<a> = NumberField(x^3 + 2)
            sage: magma(K) is magma(K)        # optional - magma
            True
            sage: magma2 = Magma()
            sage: magma(K) is magma2(K)       # optional - magma
            False
        """
        return repr(self)  # default

    def _macaulay2_(self, G=None):
        if G is None:
            import sage.interfaces.macaulay2
            G = sage.interfaces.macaulay2.macaulay2
        return self._interface_(G)

    def _macaulay2_init_(self):
        import sage.interfaces.macaulay2
        I = sage.interfaces.macaulay2.macaulay2
        return self._interface_init_(I)

    def _maple_(self, G=None):
        if G is None:
            import sage.interfaces.maple
            G = sage.interfaces.maple.maple
        return self._interface_(G)

    def _maple_init_(self):
        import sage.interfaces.maple
        I = sage.interfaces.maple.maple
        return self._interface_init_(I)

    def _mathematica_(self, G=None):
        if G is None:
            import sage.interfaces.mathematica
            G = sage.interfaces.mathematica.mathematica
        return self._interface_(G)

    def _mathematica_init_(self):
        import sage.interfaces.mathematica
        I = sage.interfaces.mathematica.mathematica
        return self._interface_init_(I)

    def _octave_(self, G=None):
        if G is None:
            import sage.interfaces.octave
            G = sage.interfaces.octave.octave
        return self._interface_(G)

    def _octave_init_(self):
        import sage.interfaces.octave
        I = sage.interfaces.octave.octave
        return self._interface_init_(I)

    def _r_init_(self):
        """
        Return default string expression that evaluates in R to this
        object.

        OUTPUT:

        - string

        EXAMPLES::

            sage: a = 2/3
            sage: a._r_init_()
            '2/3'
        """
        import sage.interfaces.r
        I = sage.interfaces.r.r
        return self._interface_init_(I)

    def _singular_(self, G=None, have_ring=False):
        if G is None:
            import sage.interfaces.singular
            G = sage.interfaces.singular.singular
        return self._interface_(G)

    def _singular_init_(self, have_ring=False):
        import sage.interfaces.singular
        I = sage.interfaces.singular.singular
        return self._interface_init_(I)

    # PARI (slightly different, since is via C library, hence instance is unique)
    def _pari_(self):
        if self._interface_is_cached_():
            try:
                return self.__pari
            except AttributeError:
                pass
        from sage.libs.pari.all import pari
        x = pari(self._pari_init_())
        if self._interface_is_cached_():
            try:
                self.__pari = x
            except AttributeError:
                # do this because C-extension class won't have a __pari attribute.
                pass
        return x

    def _pari_init_(self):
        from sage.interfaces.gp import gp
        return self._interface_init_(gp)

##################################################################

def load(*filename, compress=True, verbose=True):
    r"""
    Load Sage object from the file with name filename, which will have
    an ``.sobj`` extension added if it doesn't have one.  Or, if the input
    is a filename ending in ``.py``, ``.pyx``, ``.sage``, ``.spyx``,
    ``.f``, ``.f90`` or ``.m``, load that file into the current running
    session.

    Loaded files are not loaded into their own namespace, i.e., this is
    much more like Python's ``execfile`` than Python's ``import``.

    This function also loads a ``.sobj`` file over a network by
    specifying the full URL.  (Setting ``verbose = False`` suppresses
    the loading progress indicator.)

    Finally, if you give multiple positional input arguments, then all
    of those files are loaded, or all of the objects are loaded and a
    list of the corresponding loaded objects is returned.

    EXAMPLE::

        sage: u = 'http://sage.math.washington.edu/home/was/db/test.sobj'
        sage: s = load(u)                                                  # optional - internet
        Attempting to load remote file: http://sage.math.washington.edu/home/was/db/test.sobj
        Loading: [.]
        sage: s                                                            # optional - internet
        'hello SAGE'

    We test loading a file or multiple files or even mixing loading files and objects::

        sage: t = tmp_filename(ext='.py')
        sage: open(t,'w').write("print('hello world')")
        sage: load(t)
        hello world
        sage: load(t,t)
        hello world
        hello world
        sage: t2 = tmp_filename(); save(2/3,t2)
        sage: load(t,t,t2)
        hello world
        hello world
        [None, None, 2/3]

    Files with a ``.sage`` extension are preparsed. Also note that we
    can access global variables::

        sage: t = tmp_filename(ext=".sage")
        sage: with open(t, 'w') as f:
        ....:     f.write("a += Mod(2/3, 11)")  # This evaluates to Mod(8, 11)
        sage: a = -1
        sage: load(t)
        sage: a
        7

    We can load Fortran files::

        sage: code = '      subroutine hello\n         print *, "Hello World!"\n      end subroutine hello\n'
        sage: t = tmp_filename(ext=".F")
        sage: open(t, 'w').write(code)
        sage: load(t)
        sage: hello
        <fortran object>
    """
    import sage.repl.load
    if len(filename) != 1:
        v = [load(file, compress=compress, verbose=verbose) for file in filename]
        # Return v if one of the filenames refers to an object and not
        # a loadable filename.
        for file in filename:
            if not sage.repl.load.is_loadable_filename(file):
                return v
        return

    filename = filename[0]

    if sage.repl.load.is_loadable_filename(filename):
        sage.repl.load.load(filename, globals())
        return

    ## Check if filename starts with "http://" or "https://"
    lower = filename.lower()
    if lower.startswith("http://") or lower.startswith("https://"):
        from sage.misc.remote_file import get_remote_file
        filename = get_remote_file(filename, verbose=verbose)
        tmpfile_flag = True
    else:
        tmpfile_flag = False
        filename = process(filename)

    ## Load file by absolute filename
    X = loads(open(filename).read(), compress=compress)
    try:
        X._default_filename = os.path.abspath(filename)
    except AttributeError:
        pass

    ## Delete the tempfile, if it exists
    if tmpfile_flag == True:
        os.unlink(filename)

    return X


def save(obj, filename=None, compress=True, **kwds):
    """
    Save ``obj`` to the file with name ``filename``, which will have an
    ``.sobj`` extension added if it doesn't have one and if ``obj``
    doesn't have its own ``save()`` method, like e.g. Python tuples.

    For image objects and the like (which have their own ``save()`` method),
    you may have to specify a specific extension, e.g. ``.png``, if you
    don't want the object to be saved as a Sage object (or likewise, if
    ``filename`` could be interpreted as already having some extension).

    .. WARNING::

       This will *replace* the contents of the file if it already exists.

    EXAMPLES::

        sage: a = matrix(2, [1,2,3,-5/2])
        sage: objfile = os.path.join(SAGE_TMP, 'test.sobj')
        sage: objfile_short = os.path.join(SAGE_TMP, 'test')
        sage: save(a, objfile)
        sage: load(objfile_short)
        [   1    2]
        [   3 -5/2]
        sage: E = EllipticCurve([-1,0])
        sage: P = plot(E)
        sage: save(P, objfile_short)   # saves the plot to "test.sobj"
        sage: save(P, filename=os.path.join(SAGE_TMP, "sage.png"), xmin=-2)
        sage: save(P, os.path.join(SAGE_TMP, "filename.with.some.wrong.ext"))
        Traceback (most recent call last):
        ...
        ValueError: allowed file extensions for images are '.eps', '.pdf', '.pgf', '.png', '.ps', '.sobj', '.svg'!
        sage: print(load(objfile))
        Graphics object consisting of 2 graphics primitives
        sage: save("A python string", os.path.join(SAGE_TMP, 'test'))
        sage: load(objfile)
        'A python string'
        sage: load(objfile_short)
        'A python string'

    TESTS:

    Check that :trac:`11577` is fixed::

        sage: filename = os.path.join(SAGE_TMP, "foo.bar")  # filename containing a dot
        sage: save((1,1),filename)            # saves tuple to "foo.bar.sobj"
        sage: load(filename)
        (1, 1)
    """
    # Add '.sobj' if the filename currently has no extension
    # and also if the object doesn't have its own save() method (#11577)
    # since we otherwise assume below it is an image object:
    if not os.path.splitext(filename)[1] or not hasattr(obj, "save"):
        filename = process(filename)

    if filename.endswith('.sobj'):
        try:
            obj.save(filename=filename, compress=compress, **kwds)
        except (AttributeError, RuntimeError, TypeError):
            s = cPickle.dumps(obj, protocol=2)
            if compress:
                s = comp.compress(s)
            open(process(filename), 'wb').write(s)
    else:
        # Saving an object to an image file.
        obj.save(filename, **kwds)

def dumps(obj, compress=True):
    """
    Dump obj to a string s.  To recover obj, use ``loads(s)``.

    .. SEEALSO:: :func:`dumps`

    EXAMPLES::

        sage: a = 2/3
        sage: s = dumps(a)
        sage: len(s)
        49
        sage: loads(s)
        2/3
    """
    if make_pickle_jar:
        picklejar(obj)
    try:
        return obj.dumps(compress)
    except (AttributeError, RuntimeError, TypeError):
        if compress:
            return comp.compress(cPickle.dumps(obj, protocol=2))
        else:
            return cPickle.dumps(obj, protocol=2)

# This is used below, and also by explain_pickle.py
unpickle_override = {}

def register_unpickle_override(module, name, callable, call_name=None):
    r"""
    Python pickles include the module and class name of classes.
    This means that rearranging the Sage source can invalidate old
    pickles.  To keep the old pickles working, you can call
    register_unpickle_override with an old module name and class name,
    and the Python callable (function, class with __call__ method, etc.)
    to use for unpickling.  (If this callable is a value in some module,
    you can specify the module name and class name, for the benefit of
    :func:`~sage.misc.explain_pickle.explain_pickle` when called with ``in_current_sage=True``).)

    EXAMPLES:

    Imagine that there used to be an ``old_integer`` module and old
    pickles essentially trying to do the following::

        sage: unpickle_global('sage.rings.old_integer', 'OldInteger')
        Traceback (most recent call last):
        ...
        ImportError: cannot import OldInteger from sage.rings.old_integer, call register_unpickle_override('sage.rings.old_integer', 'OldInteger', ...) to fix this

    After following the advice from the error message, unpickling
    works::

        sage: from sage.structure.sage_object import register_unpickle_override
        sage: register_unpickle_override('sage.rings.old_integer', 'OldInteger', Integer)
        sage: unpickle_global('sage.rings.old_integer', 'OldInteger')
        <type 'sage.rings.integer.Integer'>

    In many cases, unpickling problems for old pickles can be resolved with a
    simple call to ``register_unpickle_override``, as in the example above and
    in many of the ``sage`` source files.  However, if the underlying data
    structure has changed significantly then unpickling may fail and it
    will be necessary to explicitly implement unpickling methods for the
    associated objects. The python pickle protocol is described in detail on the
    web and, in particular, in the `python pickling documentation`_. For example, the
    following excerpt from this documentation shows that the unpickling of
    classes is controlled by their :meth:`__setstate__` method.

    ::

        object.__setstate__(state)

            Upon unpickling, if the class also defines the method :meth:`__setstate__`, it is
            called with the unpickled state. If there is no :meth:`__setstate__` method,
            the pickled state must be a dictionary and its items are assigned to the new
            instance's dictionary. If a class defines both :meth:`getstate__` and
            :meth:`__setstate__`, the state object needn't be a dictionary and these methods
            can do what they want.

    .. _python pickling documentation: http://docs.python.org/library/pickle.html#pickle-protocol

    By implementing a :meth:`__setstate__` method for a class it should be
    possible to fix any unpickling problems for the class. As an example of what
    needs to be done, we show how to unpickle a :class:`CombinatorialObject`
    object using a class which also inherits from
    :class:`~sage.structure.element.Element`. This exact problem often arises
    when refactoring old code into the element framework. First we create a
    pickle to play with::

        sage: from sage.structure.element import Element
        sage: class SourPickle(CombinatorialObject): pass
        sage: class SweetPickle(CombinatorialObject,Element): pass
        sage: import __main__
        sage: __main__.SourPickle=SourPickle
        sage: __main__.SweetPickle=SweetPickle  # a hack to allow us to pickle command line classes
        sage: gherkin = dumps( SourPickle([1,2,3]) )

    Using :func:`register_unpickle_override` we try to sweeten our pickle, but we are unable to eat it::

        sage: from sage.structure.sage_object import register_unpickle_override
        sage: register_unpickle_override('__main__','SourPickle',SweetPickle)
        sage: loads( gherkin )
        Traceback (most recent call last):
        ...
        KeyError: 0

    The problem is that the ``SweetPickle`` has inherited a
    :meth:`~sage.structure.element.Element.__setstate__` method from
    :class:`~sage.structure.element.Element` which is not compatible with
    unpickling for :class:`CombinatorialObject`. We can fix this by explicitly
    defining a new :meth:`__setstate__` method::

        sage: class SweeterPickle(CombinatorialObject,Element):
        ....:     def __setstate__(self, state):
        ....:         if isinstance(state, dict):       # a pickle from CombinatorialObject is just its instance dictionary
        ....:             self._set_parent(Tableaux())  # this is a fudge: we need an appropriate parent here
        ....:             self.__dict__ = state
        ....:         else:
        ....:             self._set_parent(state[0])
        ....:             self.__dict__ = state[1]
        sage: __main__.SweeterPickle = SweeterPickle
        sage: register_unpickle_override('__main__','SourPickle',SweeterPickle)
        sage: loads( gherkin )
        [1, 2, 3]
        sage: loads(dumps( SweeterPickle([1,2,3]) ))   # check that pickles work for SweeterPickle
        [1, 2, 3]

    The ``state`` passed to :meth:`__setstate__` will usually be something like
    the instance dictionary of the pickled object, however, with some older
    classes such as :class:`CombinatorialObject` it will be a tuple. In general,
    the ``state`` can be any python object.  ``Sage`` provides a special tool,
    :func:`~sage.misc.explain_pickle.explain_pickle`, which can help in figuring
    out the contents of an old pickle. Here is a second example.

    ::

        sage: class A(object):
        ....:    def __init__(self,value):
        ....:        self.original_attribute = value
        ....:    def __repr__(self):
        ....:        return 'A(%s)'%self.original_attribute
        sage: class B(object):
        ....:    def __init__(self,value):
        ....:        self.new_attribute = value
        ....:    def __setstate__(self,state):
        ....:        try:
        ....:            self.new_attribute = state['new_attribute']
        ....:        except KeyError:      # an old pickle
        ....:            self.new_attribute = state['original_attribute']
        ....:    def __repr__(self):
        ....:        return 'B(%s)'%self.new_attribute
        sage: import __main__
        sage: __main__.A=A; __main__.B=B  # a hack to allow us to pickle command line classes
        sage: A(10)
        A(10)
        sage: loads( dumps(A(10)) )
        A(10)
        sage: sage.misc.explain_pickle.explain_pickle( dumps(A(10)) )
        pg_A = unpickle_global('__main__', 'A')
        si = unpickle_newobj(pg_A, ())
        pg_make_integer = unpickle_global('sage.rings.integer', 'make_integer')
        unpickle_build(si, {'original_attribute':pg_make_integer('a')})
        si
        sage: from sage.structure.sage_object import register_unpickle_override
        sage: register_unpickle_override('__main__', 'A', B)
        sage: loads( dumps(A(10)) )
        B(10)
        sage: loads( dumps(B(10)) )
        B(10)

    Pickling for python classes and extension classes, such as cython, is
    different -- again this is discussed in the `python pickling documentation`_. For the
    unpickling of extension classes you need to write a :meth:`__reduce__`
    method which typically returns a tuple ``(f, args,...)`` such that
    ``f(*args)`` returns (a copy of) the original object. The following code
    snippet is the :meth:`~sage.rings.integer.Integer.__reduce__` method from
    :class:`sage.rings.integer.Integer`.

    .. code-block:: cython

        def __reduce__(self):
            'Including the documentation properly causes a doc-test failure so we include it as a comment:'
            #* '''
            #* This is used when pickling integers.
            #*
            #* EXAMPLES::
            #*
            #*     sage: n = 5
            #*     sage: t = n.__reduce__(); t
            #*     (<built-in function make_integer>, ('5',))
            #*     sage: t[0](*t[1])
            #*     5
            #*     sage: loads(dumps(n)) == n
            #*     True
            #* '''
            # This single line below took me HOURS to figure out.
            # It is the *trick* needed to pickle Cython extension types.
            # The trick is that you must put a pure Python function
            # as the first argument, and that function must return
            # the result of unpickling with the argument in the second
            # tuple as input. All kinds of problems happen
            # if we don't do this.
            return sage.rings.integer.make_integer, (self.str(32),)

    """
    unpickle_override[(module,name)] = (callable, call_name)

def unpickle_global(module, name):
    r"""
    Given a module name and a name within that module (typically a class
    name), retrieve the corresponding object.  This normally just looks
    up the name in the module, but it can be overridden by
    register_unpickle_override.  This is used in the Sage unpickling
    mechanism, so if the Sage source code organization changes,
    register_unpickle_override can allow old pickles to continue to work.

    EXAMPLES::

        sage: from sage.structure.sage_object import unpickle_override, register_unpickle_override
        sage: unpickle_global('sage.rings.integer', 'Integer')
        <type 'sage.rings.integer.Integer'>

    Now we horribly break the pickling system::

        sage: register_unpickle_override('sage.rings.integer', 'Integer', Rational, call_name=('sage.rings.rational', 'Rational'))
        sage: unpickle_global('sage.rings.integer', 'Integer')
        <type 'sage.rings.rational.Rational'>

    and we reach into the internals and put it back::

        sage: del unpickle_override[('sage.rings.integer', 'Integer')]
        sage: unpickle_global('sage.rings.integer', 'Integer')
        <type 'sage.rings.integer.Integer'>
        
    A meaningful error message with resolution instructions is displayed for
    old pickles that accidentally got broken because a class or entire module
    was moved or renamed::
    
        sage: unpickle_global('sage.all', 'some_old_class')
        Traceback (most recent call last):
        ...
        ImportError: cannot import some_old_class from sage.all, call
        register_unpickle_override('sage.all', 'some_old_class', ...)
        to fix this
        
        sage: unpickle_global('sage.some_old_module', 'some_old_class')
        Traceback (most recent call last):
        ...
        ImportError: cannot import some_old_class from sage.some_old_module, call
        register_unpickle_override('sage.some_old_module', 'some_old_class', ...)
        to fix this
    """
    unpickler = unpickle_override.get((module, name))
    if unpickler is not None:
        return unpickler[0]

    def error():
        raise ImportError("cannot import {1} from {0}, call "
            "register_unpickle_override({0!r}, {1!r}, ...) to fix this".format(
                module, name))

    mod = sys_modules.get(module)
    if mod is not None:
        try:
            return getattr(mod, name)
        except AttributeError:
            error()
    try:
        __import__(module)
    except ImportError:
        error()
    mod = sys_modules[module]
    return getattr(mod, name)

def loads(s, compress=True):
    """
    Recover an object x that has been dumped to a string s
    using ``s = dumps(x)``.

    .. SEEALSO:: :func:`dumps`

    EXAMPLES::

        sage: a = matrix(2, [1,2,3,-4/3])
        sage: s = dumps(a)
        sage: loads(s)
        [   1    2]
        [   3 -4/3]

    If compress is True (the default), it will try to decompress
    the data with zlib and with bz2 (in turn); if neither succeeds,
    it will assume the data is actually uncompressed.  If compress=False
    is explicitly specified, then no decompression is attempted.

    ::

        sage: v = [1..10]
        sage: loads(dumps(v, compress=False)) == v
        True
        sage: loads(dumps(v, compress=False), compress=True) == v
        True
        sage: loads(dumps(v, compress=True), compress=False)
        Traceback (most recent call last):
        ...
        UnpicklingError: invalid load key, 'x'.
    """
    if not isinstance(s, str):
        raise TypeError("s must be a string")
    if compress:
        try:
            s = comp.decompress(s)
        except Exception as msg1:
            try:
                s = comp_other.decompress(s)
            except Exception as msg2:
                # Maybe data is uncompressed?
                pass

    unpickler = cPickle.Unpickler(StringIO(s))
    unpickler.find_global = unpickle_global

    return unpickler.load()


cdef bint make_pickle_jar = 'SAGE_PICKLE_JAR' in os.environ

def picklejar(obj, dir=None):
    """
    Create pickled sobj of obj in dir, with name the absolute value of
    the hash of the pickle of obj.  This is used in conjunction with
    :func:`unpickle_all`.

    To use this to test the whole Sage library right now, set the
    environment variable ``SAGE_PICKLE_JAR``, which will make it so dumps
    will by default call picklejar with the default dir.  Once you do
    that and doctest Sage, you'll find that the ``SAGE_ROOT``/tmp/
    contains a bunch of pickled objects along with corresponding txt
    descriptions of them.  Use the :func:`unpickle_all` to see if they unpickle
    later.

    INPUT:

    - ``obj`` -- a pickleable object

    - ``dir`` -- a string or None; if None then ``dir`` defaults to
      ``SAGE_ROOT/tmp/pickle_jar``

    EXAMPLES::

        sage: dir = tmp_dir()
        sage: sage.structure.sage_object.picklejar(1, dir)
        sage: sage.structure.sage_object.picklejar('test', dir)
        sage: len(os.listdir(dir))   # Two entries (sobj and txt) for each object
        4

    TESTS:

    Test an unaccessible directory::

        sage: import os
        sage: os.chmod(dir, 0o000)
        sage: try:
        ....:     uid = os.getuid()
        ....: except AttributeError:
        ....:     uid = -1
        sage: if uid==0:
        ....:     raise OSError('You must not run the doctests as root, geez!')
        ....: else:
        ....:     sage.structure.sage_object.picklejar(1, dir + '/noaccess')
        Traceback (most recent call last):
        ...
        OSError: ...
        sage: os.chmod(dir, 0o755)
    """
    if dir is None:
        dir = os.environ['SAGE_ROOT'] + '/tmp/pickle_jar/'
    try:
        os.makedirs(dir)
    except OSError as err:
        # It is not an error if the directory exists
        import errno
        if not err.errno == errno.EEXIST:
            raise

    s = comp.compress(cPickle.dumps(obj,protocol=2))

    typ = str(type(obj))
    name = ''.join([x if (x.isalnum() or x == '_') else '_' for x in typ])
    base = '%s/%s'%(dir, name)
    if os.path.exists(base):
        i = 0
        while os.path.exists(base + '-%s'%i):
            i += 1
        base += '-%s'%i

    open(base + '.sobj', 'wb').write(s)
    txt = "type(obj) = %s\n"%typ
    import sage.version
    txt += "version = %s\n"%sage.version.version
    txt += "obj =\n'%s'\n"%str(obj)

    open(base + '.txt', 'w').write(txt)

def unpickle_all(dir = None, debug=False, run_test_suite=False):
    """
    Unpickle all sobj's in the given directory, reporting failures as
    they occur.  Also printed the number of successes and failure.

    INPUT:

    - ``dir`` -- a string; the name of a directory (or of a .tar.bz2
      file that decompresses to a directory) full of pickles.
      (default: the standard pickle jar)
    - ``debug`` -- a boolean (default: False)
      whether to report a stacktrace in case of failure
    - ``run_test_suite`` -- a boolean (default: False)
      whether to run ``TestSuite(x).run()`` on the unpickled objects

    EXAMPLES::

        sage: dir = tmp_dir()
        sage: sage.structure.sage_object.picklejar('hello', dir)
        sage: sage.structure.sage_object.unpickle_all(dir)
        Successfully unpickled 1 objects.
        Failed to unpickle 0 objects.

    When run with no arguments :meth:`unpickle_all` tests that all of the
    "standard" pickles stored in the pickle_jar at
    ``SAGE_SHARE/sage/ext/pickle_jar/pickle_jar.tar.bz2`` can be unpickled.

    ::

        sage: sage.structure.sage_object.unpickle_all()  # (4s on sage.math, 2011)
        doctest:... DeprecationWarning: ...
        See http://trac.sagemath.org/... for details.
        Successfully unpickled ... objects.
        Failed to unpickle 0 objects.

    Check that unpickling a second time works (see :trac:`5838`)::

        sage: sage.structure.sage_object.unpickle_all()
        Successfully unpickled ... objects.
        Failed to unpickle 0 objects.

    When it is not possible to unpickle a pickle in the pickle_jar then
    :meth:`unpickle_all` prints the following error message which warns against removing
    pickles from the pickle_jar and directs the user towards
    :meth:`register_unpickle_override`. The following code intentionally
    breaks a pickle to illustrate this::

        sage: from sage.structure.sage_object import register_unpickle_override, unpickle_all, unpickle_global
        sage: class A(CombinatorialObject,sage.structure.element.Element):
        ....:     pass  # to break a pickle
        sage: tableau_unpickler=unpickle_global('sage.combinat.tableau','Tableau_class')
        sage: register_unpickle_override('sage.combinat.tableau','Tableau_class',A) # breaking the pickle
        sage: unpickle_all()  # todo: not tested
        ...
        Failed:
        _class__sage_combinat_crystals_affine_AffineCrystalFromClassicalAndPromotion_with_category_element_class__.sobj
        _class__sage_combinat_crystals_tensor_product_CrystalOfTableaux_with_category_element_class__.sobj
        _class__sage_combinat_crystals_tensor_product_TensorProductOfCrystalsWithGenerators_with_category__.sobj
        _class__sage_combinat_tableau_Tableau_class__.sobj
        ----------------------------------------------------------------------
        ** This error is probably due to an old pickle failing to unpickle.
        ** See sage.structure.sage_object.register_unpickle_override for
        ** how to override the default unpickling methods for (old) pickles.
        ** NOTE: pickles should never be removed from the pickle_jar!
        ----------------------------------------------------------------------
        Successfully unpickled 583 objects.
        Failed to unpickle 4 objects.
        sage: register_unpickle_override('sage.combinat.tableau','Tableau_class',tableau_unpickler) # restore to default

    .. TODO::

        Create a custom-made ``SourPickle`` for the last example.

    If you want to find *lots* of little issues in Sage then try the following::

        sage: sage.structure.sage_object.unpickle_all(run_test_suite = True) # todo: not tested

    This runs :class:`TestSuite` tests on all objects in the Sage pickle
    jar. Some of those objects seem to unpickle properly, but do not
    pass the tests because their internal data structure is messed
    up. In most cases though it is just that their source file misses
    a TestSuite call, and therefore some misfeatures went unnoticed
    (typically Parents not implementing the ``an_element`` method).

    .. NOTE::

        Every so often the standard pickle jar should be updated by running the
        doctest suite with the environment variable ``SAGE_PICKLE_JAR`` set, then
        copying the files from ``SAGE_ROOT/tmp/pickle_jar*`` into the standard pickle
        jar.

    .. WARNING::

        Sage's pickle jar helps to ensure backward compatibility in sage. Pickles should
        **only** be removed from the pickle jar after the corresponding objects
        have been properly deprecated. Any proposal to remove pickles from the
        pickle jar should first be discussed on sage-devel.
    """
    if dir is None:
        pickle_jar=True
        from sage.env import SAGE_EXTCODE
        dir = os.path.join(SAGE_EXTCODE,'pickle_jar','pickle_jar.tar.bz2')
    else:
        pickle_jar=False
    i = 0
    j = 0
    failed = []
    tracebacks = []
    # This could use instead Python's tarfile module
    if dir.endswith('.tar.bz2'):
        # create a temporary directory
        from sage.misc.all import tmp_dir
        T = tmp_dir()
        # extract tarball to it
        os.system('cd "%s"; bunzip2 -c "%s" | tar fx - '%(T, os.path.abspath(dir)))
        # Now use the directory in the tarball instead of dir
        dir = T + "/" + os.listdir(T)[0]

    for A in sorted(os.listdir(dir)):
        if A.endswith('.sobj'):
            try:
                object = load(os.path.join(dir,A))
                if run_test_suite:
                    TestSuite(object).run(catch = False)
                i += 1
            except Exception:
                j += 1
                if run_test_suite:
                    print(" * unpickle failure: TestSuite(load('%s')).run()" % os.path.join(dir, A))
                else:
                    print(" * unpickle failure: load('%s')" % os.path.join(dir, A))
                from traceback import print_exc
                print_exc()
                failed.append(A)
                if debug:
                    tracebacks.append(sys.exc_info())

    if len(failed) > 0:
        print("Failed:\n%s" % ('\n'.join(failed)))
        if pickle_jar:
            # if we are checking the pickle_jar then print a message to alert the
            # user about what to do to fix unpickling errors
            print('-' * 70)
            print("""** This error is probably due to an old pickle failing to unpickle.
** See sage.structure.sage_object.register_unpickle_override for
** how to override the default unpickling methods for (old) pickles.
** NOTE: pickles should never be removed from the pickle_jar! """)
            print('-' * 70)
    print("Successfully unpickled %s objects." % i)
    print("Failed to unpickle %s objects." % j)
    if debug:
        return tracebacks


def make_None(*args, **kwds):
    """
    Do nothing and return ``None``. Used for overriding pickles when
    that pickle is no longer needed.

    EXAMPLES::

        sage: from sage.structure.sage_object import make_None
        sage: print(make_None(42, pi, foo='bar'))
        None
    """
    return None


# Generators is no longer used (#21382)
register_unpickle_override('sage.structure.generators', 'make_list_gens', make_None)
