# sage_setup: distribution = sagemath-gap
"""
Library Interface to GAP

This module implements a fast C library interface to GAP.
To use it, you simply call ``libgap`` (the parent of all
:class:`~sage.libs.gap.element.GapElement` instances) and use it to
convert Sage objects into GAP objects.

EXAMPLES::

    sage: a = libgap(10)
    sage: a
    10
    sage: type(a)
    <class 'sage.libs.gap.element.GapElement_Integer'>
    sage: a*a
    100
    sage: timeit('a*a')   # random output
    625 loops, best of 3: 898 ns per loop

Compared to the expect interface this is >1000 times faster::

    sage: b = gap('10')
    sage: timeit('b*b')   # random output; long time
    125 loops, best of 3: 2.05 ms per loop

If you want to evaluate GAP commands, use the :meth:`Gap.eval` method::

    sage: libgap.eval('List([1..10], i->i^2)')
    [ 1, 4, 9, 16, 25, 36, 49, 64, 81, 100 ]

not to be confused with the ``libgap`` call, which converts Sage
objects to GAP objects, for example strings to strings::

    sage: libgap('List([1..10], i->i^2)')
    "List([1..10], i->i^2)"
    sage: type(_)
    <class 'sage.libs.gap.element.GapElement_String'>

You can usually use the :meth:`~sage.libs.gap.element.GapElement.sage`
method to convert the resulting GAP element back to its Sage
equivalent::

    sage: a.sage()
    10
    sage: type(_)
    <class 'sage.rings.integer.Integer'>

    sage: libgap.eval('5/3 + 7*E(3)').sage()                                            # needs sage.rings.number_field
    7*zeta3 + 5/3

    sage: gens_of_group = libgap.AlternatingGroup(4).GeneratorsOfGroup()
    sage: generators = gens_of_group.sage()
    sage: generators   # a Sage list of Sage permutations!
    [[2, 3, 1], [1, 3, 4, 2]]
    sage: PermutationGroup(generators).cardinality()   # computed in Sage
    12
    sage: libgap.AlternatingGroup(4).Size()            # computed in GAP
    12

We can also specify which group in Sage the permutations should
consider themselves as elements of when converted to Sage::

    sage: A4 = groups.permutation.Alternating(4)
    sage: generators = gens_of_group.sage(parent=A4); generators
    [(1,2,3), (2,3,4)]
    sage: all(gen.parent() is A4 for gen in generators)
    True

So far, the following GAP data types can be directly converted to the
corresponding Sage datatype:

#. GAP booleans ``true`` / ``false`` to Sage booleans ``True`` /
   ``False``. The third GAP boolean value ``fail`` raises a
   :exc:`ValueError`.

#. GAP integers to Sage integers.

#. GAP rational numbers to Sage rational numbers.

#. GAP cyclotomic numbers to Sage cyclotomic numbers.

#. GAP permutations to Sage permutations.

#. The GAP containers ``List`` and ``rec`` are converted to Sage
   containers ``list`` and ``dict``.  Furthermore, the
   :meth:`~sage.libs.gap.element.GapElement.sage` method is applied
   recursively to the entries.

Special support is available for the GAP container classes. GAP lists
can be used as follows::

    sage: lst = libgap([1,5,7]);  lst
    [ 1, 5, 7 ]
    sage: type(lst)
    <class 'sage.libs.gap.element.GapElement_List'>
    sage: len(lst)
    3
    sage: lst[0]
    1
    sage: [ x^2 for x in lst ]
    [1, 25, 49]
    sage: type(_[0])
    <class 'sage.libs.gap.element.GapElement_Integer'>

Note that you can access the elements of GAP ``List`` objects as you
would expect from Python (with indexing starting at 0), but the
elements are still of type
:class:`~sage.libs.gap.element.GapElement`. The other GAP container
type are records, which are similar to Python dictionaries. You can
construct them directly from Python dictionaries::

    sage: libgap({'a':123, 'b':456})
    rec( a := 123, b := 456 )

Or get them as results of computations::

    sage: rec = libgap.eval('rec(a:=123, b:=456, Sym3:=SymmetricGroup(3))')
    sage: rec['Sym3']
    Sym( [ 1 .. 3 ] )
    sage: dict(rec)
    {'Sym3': Sym( [ 1 .. 3 ] ), 'a': 123, 'b': 456}

The output is a Sage dictionary whose keys are Sage strings and whose
Values are instances of :meth:`~sage.libs.gap.element.GapElement`. So,
for example, ``rec['a']`` is not a Sage integer. To recursively
convert the entries into Sage objects, you should use the
:meth:`~sage.libs.gap.element.GapElement.sage` method::

    sage: rec.sage()
    {'Sym3': NotImplementedError('cannot construct equivalent Sage object'...),
     'a': 123,
     'b': 456}

Now ``rec['a']`` is a Sage integer. We have not implemented the
conversion of the GAP symmetric group to the Sage symmetric group yet,
so you end up with a :exc:`NotImplementedError` exception object. The
exception is returned and not raised so that you can work with the
partial result.

While we don't directly support matrices yet, you can convert them to
Gap List of Lists. These lists are then easily converted into Sage
using the recursive expansion of the
:meth:`~sage.libs.gap.element.GapElement.sage` method::

    sage: M = libgap.eval('BlockMatrix([[1,1,[[1, 2],[ 3, 4]]], [1,2,[[9,10],[11,12]]], [2,2,[[5, 6],[ 7, 8]]]],2,2)')
    sage: M
    <block matrix of dimensions (2*2)x(2*2)>
    sage: M.List()   # returns a GAP List of Lists
    [ [ 1, 2, 9, 10 ], [ 3, 4, 11, 12 ], [ 0, 0, 5, 6 ], [ 0, 0, 7, 8 ] ]
    sage: M.List().sage()   # returns a Sage list of lists
    [[1, 2, 9, 10], [3, 4, 11, 12], [0, 0, 5, 6], [0, 0, 7, 8]]
    sage: matrix(ZZ, _)
    [ 1  2  9 10]
    [ 3  4 11 12]
    [ 0  0  5  6]
    [ 0  0  7  8]


Using the GAP C library from Cython
===================================

.. TODO:: Expand the following text

   We are using the GAP API provided by the GAP project since
   GAP 4.10.

AUTHORS:

  - William Stein, Robert Miller (2009-06-23): first version
  - Volker Braun, Dmitrii Pasechnik, Ivan Andrus (2011-03-25, Sage Days 29):
    almost complete rewrite; first usable version.
  - Volker Braun (2012-08-28, GAP/Singular workshop): update to
    gap-4.5.5, make it ready for public consumption.
  - Dima Pasechnik (2018-09-18, GAP Days): started the port to native
    libgap API
"""

###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   https://www.gnu.org/licenses/
###############################################################################


##############################################################################
#
#  If you want to add support for converting some GAP datatype to its
#  Sage equivalent, you have two options. Either
#
#  1. add an if-clause to GapElement.sage(). This is the easiest
#     option and you should probably start with it first
#
#  2. Subclass GapElement to GapElement_Mydatatype. In that case you
#     need to write the derived class, a factory function to
#     instantiate it (GapElement cannot be instantiated by __init__),
#     and add an if-clause in make_GapElement. See GapElement_List,
#     for example. The advantage of this more complicated approach is
#     that you can then add extra methods to your data type. For
#     example, GapElement_List instances can access the individual
#     elements with the usual gapelement[i] syntax.
#
# TODO
#
# Not all GAP types have their own GapElement. We should wrap more
# stuff. Talk to me (Volker) if you want to work on that.
#
##############################################################################

from pathlib import Path

from sage.libs.gap.gap_includes cimport *
from sage.libs.gap.util cimport *
from sage.libs.gap.element cimport *

from sage.cpython.string cimport str_to_bytes
from sage.structure.parent cimport Parent
from sage.structure.element cimport Vector
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.misc.randstate cimport current_randstate


############################################################################
### Gap  ###################################################################
############################################################################
# The libGap interpreter object Gap is the parent of the GapElements


class Gap(Parent):
    r"""
    The libgap interpreter object.

    .. NOTE::

        This object must be instantiated exactly once by the
        libgap. Always use the provided ``libgap`` instance, and never
        instantiate :class:`Gap` manually.

    EXAMPLES::

        sage: libgap.eval('SymmetricGroup(4)')
        Sym( [ 1 .. 4 ] )

    TESTS::

        sage: TestSuite(libgap).run(skip=['_test_category', '_test_elements', '_test_pickling'])
    """

    Element = GapElement

    def _coerce_map_from_(self, S):
        """
        Whether a coercion from `S` exists.

        INPUT / OUTPUT: see :mod:`sage.structure.parent`

        EXAMPLES::

            sage: libgap.has_coerce_map_from(ZZ)
            True
            sage: libgap.has_coerce_map_from(CyclotomicField(5)['x','y'])               # needs sage.rings.number_field
            True
        """
        return True

    def _element_constructor_(self, x):
        r"""
        Construct elements of this parent class.

        INPUT:

        - ``x`` -- anything that defines a GAP object

        OUTPUT: a :class:`GapElement`

        EXAMPLES::

            sage: libgap(0)   # indirect doctest
            0
            sage: libgap(ZZ(0))
            0
            sage: libgap(int(0))
            0
            sage: libgap(vector((0,1,2)))
            [ 0, 1, 2 ]
            sage: libgap(vector((1/3,2/3,4/5)))
            [ 1/3, 2/3, 4/5 ]
            sage: libgap(vector((1/3, 0.8, 3)))
            [ 0.333333, 0.8, 3. ]
        """
        initialize()
        if isinstance(x, GapElement):
            return x
        elif isinstance(x, (list, tuple, Vector)):
            return make_GapElement_List(self, make_gap_list(x))
        elif isinstance(x, dict):
            return make_GapElement_Record(self, make_gap_record(x))
        elif isinstance(x, bool):
            # attention: must come before int
            return make_GapElement_Boolean(self, GAP_True if x else GAP_False)
        elif isinstance(x, int):
            return make_GapElement_Integer(self, make_gap_integer(x))
        elif isinstance(x, str):
            return make_GapElement_String(self, make_gap_string(x))
        elif isinstance(x, Path):
            return make_GapElement_String(self, make_gap_string(str(x)))
        else:
            try:
                return x._libgap_()
            except AttributeError:
                pass
            x = str(x._gap_init_())
            return make_any_gap_element(self, gap_eval(x))

    def _construct_matrix(self, M):
        """
        Construct a LibGAP matrix.

        INPUT:

        - ``M`` -- a matrix

        OUTPUT:

        A GAP matrix, that is, a list of lists with entries over a
        common ring.

        EXAMPLES::

            sage: M = libgap._construct_matrix(identity_matrix(ZZ,2)); M
            [ [ 1, 0 ], [ 0, 1 ] ]
            sage: M.IsMatrix()
            true

            sage: M = libgap(identity_matrix(ZZ,2)); M  # syntactic sugar
            [ [ 1, 0 ], [ 0, 1 ] ]
            sage: M.IsMatrix()
            true

            sage: M = libgap(matrix(GF(3),2,2,[4,5,6,7])); M
            [ [ Z(3)^0, Z(3) ], [ 0*Z(3), Z(3)^0 ] ]
            sage: M.IsMatrix()
            true

            sage: x = polygen(QQ, 'x')
            sage: M = libgap(matrix(QQ['x'],2,2,[x,5,6,7])); M
            [ [ x, 5 ], [ 6, 7 ] ]
            sage: M.IsMatrix()
            true

        TESTS:

        We gracefully handle the case that the conversion fails (:issue:`18039`)::

            sage: F.<a> = GF(9, modulus='first_lexicographic')                          # needs sage.rings.finite_rings
            sage: libgap(Matrix(F, [[a]]))                                              # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            NotImplementedError: conversion of (Givaro) finite field element to GAP
            not implemented except for fields defined by Conway polynomials.
        """
        ring = M.base_ring()
        try:
            gap_ring = self(ring)
        except ValueError:
            raise TypeError('base ring is not supported by GAP')
        M_list = map(list, M.rows())
        return make_GapElement_List(self, make_gap_matrix(M_list, gap_ring))

    def eval(self, gap_command):
        """
        Evaluate a gap command and wrap the result.

        INPUT:

        - ``gap_command`` -- string containing a valid gap command
          without the trailing semicolon

        OUTPUT: a :class:`GapElement`

        EXAMPLES::

            sage: libgap.eval('0')
            0
            sage: libgap.eval('"string"')
            "string"
        """
        cdef GapElement elem

        if not isinstance(gap_command, str):
            gap_command = str(gap_command._gap_init_())

        initialize()
        elem = make_any_gap_element(self, gap_eval(gap_command))

        # If the element is NULL just return None instead
        if elem.value == NULL:
            return None

        return elem

    def load_package(self, pkg):
        """
        If loading fails, raise a :exc:`RuntimeError` exception.

        TESTS::

            sage: libgap.load_package("chevie")
            Traceback (most recent call last):
            ...
            RuntimeError: Error loading GAP package chevie. You may want to
            install gap_packages SPKG.
        """
        load_package = self.function_factory('LoadPackage')
        # Note: For some reason the default package loading error messages are
        # controlled with InfoWarning and not InfoPackageLoading
        prev_infolevel = libgap.InfoLevel(libgap.InfoWarning)
        libgap.SetInfoLevel(libgap.InfoWarning, 0)
        ret = load_package(pkg)
        libgap.SetInfoLevel(libgap.InfoWarning, prev_infolevel)
        if str(ret) == 'fail':
            raise RuntimeError(f"Error loading GAP package {pkg}.  "
                               f"You may want to install gap_packages SPKG.")
        return ret

    @cached_method
    def function_factory(self, function_name):
        """
        Return a GAP function wrapper.

        This is almost the same as calling
        ``libgap.eval(function_name)``, but faster and makes it
        obvious in your code that you are wrapping a function.

        INPUT:

        - ``function_name`` -- string; the name of a GAP function

        OUTPUT:

        A function wrapper
        :class:`~sage.libs.gap.element.GapElement_Function` for the
        GAP function. Calling it from Sage is equivalent to calling
        the wrapped function from GAP.

        EXAMPLES::

            sage: libgap.function_factory('Print')
            <Gap function "Print">
        """
        initialize()
        return make_GapElement_Function(self, gap_eval(function_name))

    def set_global(self, variable, value):
        """
        Set a GAP global variable.

        INPUT:

        - ``variable`` -- string; the variable name

        - ``value`` -- anything that defines a GAP object

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: libgap.get_global('FooBar')
            1
            sage: libgap.unset_global('FooBar')
            sage: libgap.get_global('FooBar')
            NULL
        """
        is_bound = self.function_factory('IsBoundGlobal')
        bind_global = self.function_factory('BindGlobal')
        if is_bound(variable):
            self.unset_global(variable)
        bind_global(variable, value)

    def unset_global(self, variable):
        """
        Remove a GAP global variable.

        INPUT:

        - ``variable`` -- string; the variable name

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: libgap.get_global('FooBar')
            1
            sage: libgap.unset_global('FooBar')
            sage: libgap.get_global('FooBar')
            NULL
        """
        is_readonlyglobal = self.function_factory('IsReadOnlyGlobal')
        make_readwrite = self.function_factory('MakeReadWriteGlobal')
        unbind_global = self.function_factory('UnbindGlobal')
        if is_readonlyglobal(variable):
            make_readwrite(variable)
        unbind_global(variable)

    def get_global(self, variable):
        """
        Get a GAP global variable.

        INPUT:

        - ``variable`` -- string; the variable name

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement` wrapping the GAP
        output. A :exc:`ValueError` is raised if there is no such
        variable in GAP.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: libgap.get_global('FooBar')
            1
            sage: libgap.unset_global('FooBar')
            sage: libgap.get_global('FooBar')
            NULL
        """
        return make_any_gap_element(self, GAP_ValueGlobalVariable(str_to_bytes(variable)))

    def global_context(self, variable, value):
        """
        Temporarily change a global variable.

        INPUT:

        - ``variable`` -- string; the variable name

        - ``value`` -- anything that defines a GAP object

        OUTPUT: a context manager that sets/reverts the given global variable

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: with libgap.global_context('FooBar', 2):
            ....:     print(libgap.get_global('FooBar'))
            2
            sage: libgap.get_global('FooBar')
            1
        """
        from sage.libs.gap.context_managers import GlobalVariableContext
        initialize()
        return GlobalVariableContext(variable, value)

    def set_seed(self, seed=None):
        """
        Reseed the standard GAP pseudo-random sources with the given seed.

        Uses a random seed given by ``current_randstate().ZZ_seed()`` if
        ``seed=None``.  Otherwise the seed should be an integer.

        EXAMPLES::

            sage: libgap.set_seed(0)
            0
            sage: [libgap.Random(1, 10) for i in range(5)]
            [2, 3, 3, 4, 2]
        """
        if seed is None:
            seed = current_randstate().ZZ_seed()

        Reset = self.function_factory("Reset")
        Reset(self.GlobalMersenneTwister, seed)
        Reset(self.GlobalRandomSource, seed)
        return seed

    def _an_element_(self):
        r"""
        Return a :class:`GapElement`.

        OUTPUT: a :class:`GapElement`

        EXAMPLES::

            sage: libgap.an_element()   # indirect doctest
            0
        """
        return self(0)

    def zero(self):
        """
        Return (integer) zero in GAP.

        OUTPUT: a :class:`GapElement`

        EXAMPLES::

            sage: libgap.zero()
            0
        """
        return self(0)

    def one(self):
        r"""
        Return (integer) one in GAP.

        EXAMPLES::

            sage: libgap.one()
            1
            sage: parent(_)
            C library interface to GAP
        """
        return self(1)

    def __init__(self):
        r"""
        The Python constructor.

        EXAMPLES::

            sage: type(libgap)
            <class 'sage.misc.lazy_import.LazyImport'>
            sage: type(libgap._get_object())
            <class 'sage.libs.gap.libgap.Gap'>
        """
        Parent.__init__(self, base=ZZ)

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        EXAMPLES::

            sage: libgap
            C library interface to GAP
        """
        return 'C library interface to GAP'

    @cached_method
    def __dir__(self):
        """
        Customize tab completion.

        EXAMPLES::

           sage: 'OctaveAlgebra' in dir(libgap)
           True
        """
        from sage.libs.gap.gap_globals import common_gap_globals
        return dir(self.__class__) + sorted(common_gap_globals)

    def __getattr__(self, name):
        r"""
        The attributes of the Gap object are the Gap functions, and in some
        cases other global variables from GAP.

        INPUT:

        - ``name`` -- string; the name of the GAP function you want to
          call

        OUTPUT:

        A :class:`GapElement`. A :exc:`AttributeError` is raised
        if there is no such function or global variable.

        EXAMPLES::

            sage: libgap.List
            <Gap function "List">
            sage: libgap.GlobalRandomSource
            <RandomSource in IsGAPRandomSource>
        """
        if name in dir(self.__class__):
            return getattr(self.__class__, name)

        try:
            g = self.eval(name)
        except ValueError:
            raise AttributeError(f'No such attribute: {name}.')

        self.__dict__[name] = g
        return g

    def show(self):
        """
        Return statistics about the GAP owned object list.

        This includes the total memory allocated by GAP as returned by
        ``libgap.eval('TotalMemoryAllocated()'), as well as garbage collection
        / object count statistics as returned by
        ``libgap.eval('GasmanStatistics')``, and finally the total number of
        GAP objects held by Sage as :class:`~sage.libs.gap.element.GapElement`
        instances.

        The value ``livekb + deadkb`` will roughly equal the total memory
        allocated for GAP objects (see
        ``libgap.eval('TotalMemoryAllocated()')``).

        .. NOTE::

            Slight complication is that we want to do it without accessing
            libgap objects, so we don't create new GapElements as a side
            effect.

        EXAMPLES::

            sage: a = libgap(123)
            sage: b = libgap(456)
            sage: c = libgap(789)
            sage: del b
            sage: libgap.collect()
            sage: libgap.show()  # random output
            {'gasman_stats': {'full': {'cumulative': 110,
               'deadbags': 321400,
               'deadkb': 12967,
               'freekb': 15492,
               'livebags': 396645,
               'livekb': 37730,
               'time': 110,
               'totalkb': 65536},
              'nfull': 1,
              'npartial': 1},
             'nelements': 23123,
             'total_alloc': 3234234}
        """
        d = {'nelements': self.count_GAP_objects()}
        d['total_alloc'] = self.eval('TotalMemoryAllocated()').sage()
        d['gasman_stats'] = self.eval('GasmanStatistics()').sage()
        return d

    def count_GAP_objects(self):
        """
        Return the number of GAP objects that are being tracked by
        GAP.

        OUTPUT: integer

        EXAMPLES::

            sage: libgap.count_GAP_objects()   # random output
            5
        """
        return len(get_owned_objects())

    def collect(self):
        """
        Manually run the garbage collector.

        EXAMPLES::

            sage: a = libgap(123)
            sage: del a
            sage: libgap.collect()
        """
        initialize()
        GAP_CollectBags(1)


libgap = Gap()
