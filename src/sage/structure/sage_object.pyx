# sage_setup: distribution = sagemath-objects
r"""
Abstract base class for Sage objects
"""

from sage.misc.persist import (_base_dumps, _base_save,
                               register_unpickle_override, make_None)

from sage.misc.lazy_import import LazyImport

# NOTE: These imports are just for backwards-compatibility
loads = LazyImport('sage.misc.persist', 'loads', deprecation=25153)
dumps = LazyImport('sage.misc.persist', 'dumps', deprecation=25153)
save = LazyImport('sage.misc.persist', 'save', deprecation=25153)
load = LazyImport('sage.misc.persist', 'load', deprecation=25153)
unpickle_all = LazyImport('sage.misc.persist', 'unpickle_all',
                          deprecation=25153)
unpickle_global = LazyImport('sage.misc.persist', 'unpickle_global',
                             deprecation=25153)
unpickle_override = LazyImport('sage.misc.persist', 'unpickle_override',
                               deprecation=25153)


# Generators is no longer used (#21382)
register_unpickle_override('sage.structure.generators', 'make_list_gens',
                           make_None)


__all__ = ['SageObject']


# The _interface_init_ for these interfaces takes the interface as argument
_interface_init_with_interface = set(['magma', 'macaulay2'])


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
        Change ``self`` so it prints as x, where x is a string.

        If x is ``None``, the existing custom name is removed.

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
            sage: g.rename(None)
            sage: g
            x^3 + x - 5

        Real numbers are not Python classes, so rename is not supported::

            sage: a = 3.14
            sage: type(a)                                                               # needs sage.rings.real_mpfr
            <... 'sage.rings.real_mpfr.RealLiteral'>
            sage: a.rename('pi')                                                        # needs sage.rings.real_mpfr
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support renaming: 3.14000000000000

        .. NOTE::

           The reason C-extension types are not supported by default
           is if they were then every single one would have to carry
           around an extra attribute, which would be slower and waste
           a lot of memory.

           To support them for a specific class, add a
           ``cdef public _SageObject__custom_name`` attribute.
        """
        if x is None:
            self.reset_name()
        else:
            try:
                # TODO: after dropping support for Cython < 3.0.0, all
                # the self._SageObject__custom_name in this class can be
                # changed to self.__custom_name
                self._SageObject__custom_name = str(x)
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
        if hasattr(self, '_SageObject__custom_name'):
            del self._SageObject__custom_name

    def get_custom_name(self):
        """
        Return the custom name of this object, or ``None`` if it is not
        renamed.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: P.get_custom_name() is None
            True
            sage: P.rename('A polynomial ring')
            sage: P.get_custom_name()
            'A polynomial ring'
            sage: P.reset_name()
            sage: P.get_custom_name() is None
            True
        """
        try:
            return self._SageObject__custom_name
        except AttributeError:
            return None

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

        If there is no ``_repr_`` method defined, we fall back to the
        super class (typically ``object``)::

            sage: from sage.structure.sage_object import SageObject
            sage: S = SageObject()
            sage: S
            <sage.structure.sage_object.SageObject object at ...>
        """
        try:
            name = self._SageObject__custom_name
            if name is not None:
                return name
        except AttributeError:
            pass
        try:
            reprfunc = self._repr_
        except AttributeError:
            return super().__repr__()
        return reprfunc()

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

            sage: result = ascii_art(integral(exp(x+x^2)/(x+1), x))                     # needs sage.symbolic
            ...
            sage: result                                                                # needs sage.symbolic
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

            sage: # needs sage.combinat
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

            sage: unicode_art(integral(exp(x+x^2)/(x+1), x))                            # needs sage.symbolic
            ⌠
            ⎮   2
            ⎮  x  + x
            ⎮ ℯ
            ⎮ ─────── dx
            ⎮  x + 1
            ⌡


        Alternatively, you can use the ``%display ascii_art/simple`` magic to
        switch all output to ASCII art and back::

            sage: # needs sage.combinat
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

        Check that breakpoints and baseline are preserved (:issue:`29202`)::

            sage: # needs sage.groups
            sage: F = FreeAbelianMonoid(index_set=ZZ)
            sage: f = prod(F.gen(i) for i in range(5))
            sage: s, t = ascii_art(f), unicode_art(f)
            sage: s._breakpoints == t._breakpoints and s._baseline == t._baseline
            True
        """
        from sage.typeset.unicode_art import UnicodeArt
        s = self._ascii_art_()
        lines = [unicode(z) for z in s]
        return UnicodeArt(lines, s._breakpoints, s._baseline)

    def __hash__(self):
        r"""
        Not implemented: mutable objects inherit from this class.

        EXAMPLES::

            sage: hash(SageObject())
            Traceback (most recent call last):
            ...
            TypeError: <... 'sage.structure.sage_object.SageObject'> is not hashable
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

            sage: # needs sage.rings.padics
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

            sage: # needs sage.rings.padics
            sage: hash(b)
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'sage.rings.padics.qadic_flint_CR.qAdicCappedRelativeElement'
            sage: @cached_method
            ....: def f(x): return x==a
            sage: f(b)
            True
            sage: f(c)  # if b and c were hashable, this would return True
            False
            sage: b._cache_key()
            (..., ((0, 1),), 0, 1)
            sage: c._cache_key()
            (..., ((0, 1), (1,)), 0, 20)

        An implementation must make sure that for elements ``a`` and ``b``,
        if ``a != b``, then also ``a._cache_key() != b._cache_key()``.
        In practice this means that the ``_cache_key`` should always include
        the parent as its first argument::

            sage: S.<a> = Qq(4)                                                         # needs sage.rings.padics
            sage: d = a + O(2)                                                          # needs sage.rings.padics
            sage: b._cache_key() == d._cache_key()  # this would be True if the parents were not included               # needs sage.rings.padics
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
        Save ``self`` to the given filename.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: x = SR.var("x")
            sage: f = x^3 + 5
            sage: from tempfile import NamedTemporaryFile
            sage: with NamedTemporaryFile(suffix='.sobj') as t:
            ....:     f.save(t.name)
            ....:     load(t.name)
            x^3 + 5
        """
        if filename is None:
            try:
                filename = self._default_filename
            except AttributeError:
                raise RuntimeError(
                        "no default filename, so it must be specified")

        filename = _base_save(self, filename, compress=compress)

        try:
            self._default_filename = filename
        except AttributeError:
            pass

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

            sage: from sage.misc.persist import comp
            sage: O = SageObject()
            sage: p_comp = O.dumps()
            sage: p_uncomp = O.dumps(compress=False)
            sage: comp.decompress(p_comp) == p_uncomp
            True
            sage: import pickletools
            sage: pickletools.dis(p_uncomp)
                0: \x80 PROTO      2
                2: c    GLOBAL     'sage.structure.sage_object SageObject'
               41: q    BINPUT     ...
               43: )    EMPTY_TUPLE
               44: \x81 NEWOBJ
               45: q    BINPUT     ...
               47: .    STOP
            highest protocol among opcodes = 2
        """

        return _base_dumps(self, compress=compress)

    #############################################################################
    # Category theory / structure
    #############################################################################

    def category(self):
        from sage.categories.objects import Objects
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
            AssertionError: 3 is not an instance of
            <class 'sage.categories.category.Category'>
        """
        from sage.categories.category import Category
        from sage.categories.objects import Objects
        tester = self._tester(**options)
        category = self.category()
        tester.assertIsInstance(category, Category)
        tester.assertTrue(category.is_subcategory(Objects()))
        tester.assertIn(self, category)

    def parent(self):
        """
        Return the type of ``self`` to support the coercion framework.

        EXAMPLES::

            sage: t = log(sqrt(2) - 1) + log(sqrt(2) + 1); t                            # needs sage.symbolic
            log(sqrt(2) + 1) + log(sqrt(2) - 1)
            sage: u = t.maxima_methods()                                                # needs sage.symbolic
            sage: u.parent()                                                            # needs sage.symbolic
            <class 'sage.symbolic.maxima_wrapper.MaximaWrapper'>
        """
        return type(self)

    ##########################################################################
    # Test framework
    ##########################################################################

    def _tester(self, **options):
        """
        Return a gadget attached to ``self`` providing testing utilities.

        This is used by :class:`sage.misc.sage_unittest.TestSuite` and the
        ``_test_*`` methods.

        EXAMPLES::

            sage: tester = ZZ._tester()

            sage: tester.assertTrue(1 == 1)
            sage: tester.assertTrue(1 == 0)
            Traceback (most recent call last):
            ...
            AssertionError: False is not true
            sage: tester.assertTrue(1 == 0, "this is expected to fail")
            Traceback (most recent call last):
            ...
            AssertionError:... this is expected to fail

            sage: tester.assertEqual(1, 1)
            sage: tester.assertEqual(1, 0)
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
        Check that all required methods for this object are implemented.

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

        Check that only errors triggered by ``AbstractMethod`` are caught
        (:issue:`29694`)::

            sage: class NotAbstract(SageObject):
            ....:     @lazy_attribute
            ....:     def bla(self):
            ....:         raise NotImplementedError("not implemented")
            sage: NotAbstract()._test_not_implemented_methods()
        """
        tester = self._tester(**options)
        try:
            # Disable warnings for the duration of the test
            import warnings
            warnings.filterwarnings('ignore')
            for name in dir(self):
                try:
                    getattr(self, name)
                except NotImplementedError as e:
                    if 'abstract method' in str(e):
                        tester.fail("Not implemented method: %s" % name)
                except Exception:
                    pass
        finally:
            # Restore warnings
            warnings.filters.pop(0)

    def _test_pickling(self, **options):
        """
        Check that this object can be pickled and unpickled properly.

        EXAMPLES::

            sage: ZZ._test_pickling()

        .. SEEALSO::

            :func:`dumps`, :func:`loads`

        TESTS::

            sage: class Bla(SageObject): pass
            sage: Bla()._test_pickling()
            Traceback (most recent call last):
            ...
            PicklingError: Can't pickle <class '__main__.Bla'>: attribute
            lookup ... failed

        TODO: for a stronger test, this could send the object to a
        remote Sage session, and get it back.
        """
        tester = self._tester(**options)
        from sage.misc.persist import loads, dumps
        tester.assertEqual(loads(dumps(self)), self)

    #############################################################################
    # Coercions to interface objects
    #############################################################################

    # Sage
    def _sage_(self):
        return self

    def _interface_(self, I):
        """
        Return coercion of ``self`` to an object of the interface I.

        The result of coercion is cached, unless ``self`` is a C extension
        class or ``self._interface_is_cached_()`` returns ``False``.
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
            if nm in _interface_init_with_interface:
                s = init_func(I)
            else:
                s = init_func()
        else:
            try:
                s = self._interface_init_(I)
            except Exception:
                raise NotImplementedError("coercion of object %s to %s not implemented" % (repr(self), I))
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
        Return ``True`` if the interface objects are cached.

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

    def _libgap_(self):
        from sage.libs.gap.libgap import libgap
        return libgap.eval(self)

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
        that interpreter to the Magma object corresponding to ``self``.
        This function may call the magma interpreter when it runs.

        INPUT:

        - ``magma`` -- a Magma interface

        OUTPUT: string

        EXAMPLES::

            sage: n = -3/7
            sage: n._magma_init_(magma)
            '-3/7'

        Some other examples that illustrate conversion to Magma.
        ::

            sage: # optional - magma, needs sage.symbolic
            sage: n = -3/7
            sage: m2 = Magma()
            sage: magma(n)
            -3/7
            sage: magma(n).parent()
            Magma
            sage: magma(n).parent() is m2
            False
            sage: magma(n).parent() is magma
            True

        This example illustrates caching, which happens automatically
        since K is a Python object::

            sage: # optional - magma, needs sage.symbolic
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^3 + 2)                                          # needs sage.rings.number_field
            sage: magma(K) is magma(K)
            True
            sage: magma2 = Magma()
            sage: magma(K) is magma2(K)
            False
        """
        return repr(self)  # default

    def _macaulay2_(self, G=None):
        if G is None:
            import sage.interfaces.macaulay2
            G = sage.interfaces.macaulay2.macaulay2
        return self._interface_(G)

    def _macaulay2_init_(self, macaulay2=None):
        if macaulay2 is None:
            import sage.interfaces.macaulay2
            macaulay2 = sage.interfaces.macaulay2.macaulay2
        return self._interface_init_(macaulay2)

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

    def _mathics_(self, G=None):
        if G is None:
            import sage.interfaces.mathics
            G = sage.interfaces.mathics.mathics
        return self._interface_(G)

    _mathics_init_ = _mathematica_init_

    def _octave_(self, G=None):
        if G is None:
            import sage.interfaces.octave
            G = sage.interfaces.octave.octave
        return self._interface_(G)

    def _octave_init_(self):
        import sage.interfaces.octave
        I = sage.interfaces.octave.octave
        return self._interface_init_(I)

    def _polymake_(self, G=None):
        if G is None:
            import sage.interfaces.polymake
            G = sage.interfaces.polymake.polymake
        return self._interface_(G)

    def _polymake_init_(self):
        import sage.interfaces.polymake
        I = sage.interfaces.polymake.polymake
        return self._interface_init_(I)

    def _r_init_(self):
        """
        Return default string expression that evaluates in R to this
        object.

        OUTPUT: string

        EXAMPLES::

            sage: a = 2/3
            sage: a._r_init_()                               # optional - rpy2
            '2/3'
        """
        import sage.interfaces.r
        I = sage.interfaces.r.r
        return self._interface_init_(I)

    def _singular_(self, G=None):
        if G is None:
            import sage.interfaces.singular
            G = sage.interfaces.singular.singular
        return self._interface_(G)

    def _singular_init_(self):
        import sage.interfaces.singular
        I = sage.interfaces.singular.singular
        return self._interface_init_(I)

    # PARI (slightly different, since is via C library, hence instance is unique)
    def __pari__(self):
        if self._interface_is_cached_():
            try:
                return self.__pari
            except AttributeError:
                pass
        from sage.libs.pari import pari
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
