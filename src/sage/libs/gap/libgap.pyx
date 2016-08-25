"""
libGAP shared library Interface to GAP

This module implements a fast C library interface to GAP. To use
libGAP you simply call ``libgap`` (the parent of all
:class:`~sage.libs.gap.element.GapElement` instances) and use it to
convert Sage objects into GAP objects.

EXAMPLES::

    sage: a = libgap(10)
    sage: a
    10
    sage: type(a)
    <type 'sage.libs.gap.element.GapElement_Integer'>
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
    <type 'sage.libs.gap.element.GapElement_String'>

You can usually use the :meth:`~sage.libs.gap.element.GapElement.sage`
method to convert the resulting GAP element back to its Sage
equivalent::

    sage: a.sage()
    10
    sage: type(_)
    <type 'sage.rings.integer.Integer'>

    sage: libgap.eval('5/3 + 7*E(3)').sage()
    7*zeta3 + 5/3

    sage: generators = libgap.AlternatingGroup(4).GeneratorsOfGroup().sage()
    sage: generators   # a Sage list of Sage permutations!
    [(1,2,3), (2,3,4)]
    sage: PermutationGroup(generators).cardinality()   # computed in Sage
    12
    sage: libgap.AlternatingGroup(4).Size()            # computed in GAP
    12

So far, the following GAP data types can be directly converted to the
corresponding Sage datatype:

#. GAP booleans ``true`` / ``false`` to Sage booleans ``True`` /
   ``False``. The third GAP boolean value ``fail`` raises a
   ``ValueError``.

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
    <type 'sage.libs.gap.element.GapElement_List'>
    sage: len(lst)
    3
    sage: lst[0]
    1
    sage: [ x^2 for x in lst ]
    [1, 25, 49]
    sage: type(_[0])
    <type 'sage.libs.gap.element.GapElement_Integer'>

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
    {'Sym3': NotImplementedError('cannot construct equivalent Sage object',),
     'a': 123,
     'b': 456}

Now ``rec['a']`` is a Sage integer. We have not implemented the
conversion of the GAP symmetric group to the Sage symmetric group yet,
so you end up with a ``NotImplementedError`` exception object. The
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


Using the libGAP C library from Cython
======================================

The lower-case ``libgap_foobar`` functions are ones that we added to
make the libGAP C shared library. The ``libGAP_foobar`` methods are
the original GAP methods simply prefixed with the string
``libGAP_``. The latter were originally not designed to be in a
library, so some care needs to be taken to call them.

In particular, you must call ``libgap_mark_stack_bottom()`` in every
function that calls into the libGAP C functions. The reason is that
the GAP memory manager will automatically keep objects alive that are
referenced in local (stack-allocated) variables. While convenient,
this requires to look through the stack to find anything that looks
like an address to a memory bag. But this requires vigilance against
the following pattern::

    cdef f()
      libgap_mark_stack_bottom()
      libGAP_function()

    cdef g()
      libgap_mark_stack_bottom();
      f()                #  f() changed the stack bottom marker
      libGAP_function()  #  boom

The solution is to re-order ``g()`` to first call ``f()``. In order to
catch this error, it is recommended that you wrap calls into libGAP in
``libgap_enter`` / ``libgap_exit`` blocks and not call
``libgap_mark_stack_bottom`` manually. So instead, always write

    cdef f()
      libgap_enter()
      libGAP_function()
      libgap_exit()

    cdef g()
      f()
      libgap_enter()
      libGAP_function()
      libgap_exit()

If you accidentally call ``libgap_enter()`` twice then an error
message is printed to help you debug this::

    sage: from sage.libs.gap.util import error_enter_libgap_block_twice
    sage: error_enter_libgap_block_twice()
    Traceback (most recent call last):
    ...
    RuntimeError: Entered a critical block twice

AUTHORS:

  - William Stein, Robert Miller (2009-06-23): first version
  - Volker Braun, Dmitrii Pasechnik, Ivan Andrus (2011-03-25, Sage Days 29):
    almost complete rewrite; first usable version.
  - Volker Braun (2012-08-28, GAP/Singular workshop): update to
    gap-4.5.5, make it ready for public consumption.
"""

###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
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
from __future__ import print_function

from gap_includes cimport *

from sage.structure.sage_object cimport SageObject
from sage.structure.parent cimport Parent
from sage.structure.element cimport ModuleElement, RingElement
from sage.rings.all import ZZ
from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecated_function_alias
from sage.libs.gap.element cimport *


############################################################################
### Debugging ##############################################################
############################################################################

cdef void report(libGAP_Obj bag):
    print(libGAP_TNAM_OBJ(bag), <int>libGAP_TNUM_BAG(bag), <int>libGAP_SIZE_BAG(bag))


cdef void print_gasman_objects():
    libgap_enter()
    libGAP_CallbackForAllBags(report)
    libgap_exit()




from sage.misc.lazy_import import is_during_startup
if is_during_startup():
    import sys, traceback
    print('Importing libgap during startup!')
    traceback.print_stack(None, None, sys.stdout)


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

        INPUT / OUTPUT:

        See :mod:`sage.structure.parent`.

        EXAMPLES::

            sage: libgap.has_coerce_map_from(ZZ)
            True
            sage: libgap.has_coerce_map_from(CyclotomicField(5)['x','y'])
            True
        """
        return True

    def _element_constructor_(self, x):
        r"""
        Construct elements of this parent class.

        INPUT:

        - ``x`` -- anything that defines a GAP object.

        OUTPUT:

        A :class:`GapElement`.

        EXAMPLES::

            sage: libgap(0)   # indirect doctest
            0
            sage: libgap(ZZ(0))
            0
            sage: libgap(int(0))
            0
        """
        if isinstance(x, GapElement):
            return x
        elif isinstance(x, (list, tuple)):
            return make_GapElement_List(self, make_gap_list(x))
        elif isinstance(x, dict):
            return make_GapElement_Record(self, make_gap_record(x))
        elif isinstance(x, bool):
            # attention: must come before int
            return make_GapElement_Boolean(self, libGAP_True if x else libGAP_False)
        elif isinstance(x, int):
            return make_GapElement_Integer(self, make_gap_integer(x))
        elif isinstance(x, basestring):
            return make_GapElement_String(self, make_gap_string(x))
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

        - ``M`` -- a matrix.

        OUTPUT:

        A GAP matrix, that is, a list of lists with entries over a
        common ring.

        EXAMPLES::

            sage: libgap._construct_matrix(identity_matrix(ZZ,2))
            [ [ 1, 0 ], [ 0, 1 ] ]
            sage: libgap(identity_matrix(ZZ,2))  # syntactic sugar
            [ [ 1, 0 ], [ 0, 1 ] ]
            sage: libgap(matrix(GF(3),2,2,[4,5,6,7]))
            [ [ Z(3)^0, Z(3) ], [ 0*Z(3), Z(3)^0 ] ]

        TESTS:

        We gracefully handle the case that the conversion fails (:trac:`18039`)::

            sage: F.<a> = GF(9, modulus="first_lexicographic")
            sage: libgap(Matrix(F, [[a]]))
            Traceback (most recent call last):
            ...
            NotImplementedError: conversion of (Givaro) finite field element to GAP not implemented except for fields defined by Conway polynomials.
        """
        ring = M.base_ring()
        try:
            gap_ring = self(ring)
        except ValueError:
            raise TypeError('base ring is not supported by GAP')
        M_list = map(list, M.rows())
        return make_GapElement_List(self, make_gap_list(M_list))

    def eval(self, gap_command):
        """
        Evaluate a gap command and wrap the result.

        INPUT:

        - ``gap_command`` -- a string containing a valid gap command
          without the trailing semicolon.

        OUTPUT:

        A :class:`GapElement`.

        EXAMPLES::

            sage: libgap.eval('0')
            0
            sage: libgap.eval('"string"')
            "string"
        """
        if not isinstance(gap_command, basestring):
            gap_command = str(gap_command._gap_init_())
        return make_any_gap_element(self, gap_eval(gap_command))

    @cached_method
    def function_factory(self, function_name):
        """
        Return a GAP function wrapper

        This is almost the same as calling
        ``libgap.eval(function_name)``, but faster and makes it
        obvious in your code that you are wrapping a function.

        INPUT:

        - ``function_name`` -- string. The name of a GAP function.

        OUTPUT:

        A function wrapper
        :class:`~sage.libs.gap.element.GapElement_Function` for the
        GAP function. Calling it from Sage is equivalent to calling
        the wrapped function from GAP.

        EXAMPLES::

            sage: libgap.function_factory('Print')
            <Gap function "Print">
        """
        return make_GapElement_Function(self, gap_eval(function_name))

    def set_global(self, variable, value):
        """
        Set a GAP global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        - ``value`` -- anything that defines a GAP object.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: libgap.get_global('FooBar')
            1
            sage: libgap.unset_global('FooBar')
            sage: libgap.get_global('FooBar')
            Traceback (most recent call last):
            ...
            ValueError: libGAP: Error, VAL_GVAR: No value bound to FooBar
        """
        is_bound = self.function_factory('IsBoundGlobal')
        bind_global = self.function_factory('BindGlobal')
        if is_bound(variable):
            self.unset_global(variable)
        bind_global(variable, value)

    def unset_global(self, variable):
        """
        Remove a GAP global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: libgap.get_global('FooBar')
            1
            sage: libgap.unset_global('FooBar')
            sage: libgap.get_global('FooBar')
            Traceback (most recent call last):
            ...
            ValueError: libGAP: Error, VAL_GVAR: No value bound to FooBar
        """
        make_readwrite = self.function_factory('MakeReadWriteGlobal')
        unbind_global = self.function_factory('UnbindGlobal')
        make_readwrite(variable)
        unbind_global(variable)

    def get_global(self, variable):
        """
        Get a GAP global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        OUTPUT:

        A :class:`~sage.libs.gap.element.GapElement` wrapping the GAP
        output. A ``ValueError`` is raised if there is no such
        variable in GAP.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: libgap.get_global('FooBar')
            1
            sage: libgap.unset_global('FooBar')
            sage: libgap.get_global('FooBar')
            Traceback (most recent call last):
            ...
            ValueError: libGAP: Error, VAL_GVAR: No value bound to FooBar
        """
        value_global = self.function_factory('ValueGlobal')
        return value_global(variable)

    def global_context(self, variable, value):
        """
        Temporarily change a global variable

        INPUT:

        - ``variable`` -- string. The variable name.

        - ``value`` -- anything that defines a GAP object.

        OUTPUT:

        A context manager that sets/reverts the given global variable.

        EXAMPLES::

            sage: libgap.set_global('FooBar', 1)
            sage: with libgap.global_context('FooBar', 2):
            ....:     print(libgap.get_global('FooBar'))
            2
            sage: libgap.get_global('FooBar')
            1
        """
        from sage.libs.gap.context_managers import GlobalVariableContext
        return GlobalVariableContext(variable, value)

    def _an_element_(self):
        r"""
        Return a :class:`GapElement`.

        OUTPUT:

        A :class:`GapElement`.

        EXAMPLES::

            sage: libgap.an_element()   # indirect doctest
            0
        """
        return self(0)

    def zero(self):
        """
        Return (integer) zero in GAP.

        OUTPUT:

        A :class:`GapElement`.

        EXAMPLES::

            sage: libgap.zero()
            0

        TESTS::

            sage: libgap.zero_element()
            doctest:...: DeprecationWarning: zero_element is deprecated. Please use zero instead.
            See http://trac.sagemath.org/17694 for details.
            0
        """
        return self(0)

    zero_element = deprecated_function_alias(17694, zero)

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
            <type 'sage.misc.lazy_import.LazyImport'>
            sage: type(libgap._get_object())
            <class 'sage.libs.gap.libgap.Gap'>
        """
        initialize()
        libgap_set_gasman_callback(gasman_callback)
        from sage.rings.integer_ring import ZZ
        Parent.__init__(self, base=ZZ)

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: libgap
            C library interface to GAP
        """
        return 'C library interface to GAP'

    @cached_method
    def __dir__(self):
        """
        Customize tab completion

        EXAMPLES::

           sage: 'OctaveAlgebra' in dir(libgap)
           True
        """
        from sage.libs.gap.gap_functions import common_gap_functions
        return dir(self.__class__) + list(common_gap_functions)

    def __getattr__(self, name):
        r"""
        The attributes of the Gap object are the Gap functions.

        INPUT:

        - ``name`` -- string. The name of the GAP function you want to
          call.

        OUTPUT:

        A :class:`GapElement_Function`. A ``AttributeError`` is raised
        if there is no such function.

        EXAMPLES::

            sage: libgap.List
            <Gap function "List">
        """
        if name in dir(self.__class__):
            return getattr(self.__class__, name)
        from sage.libs.gap.gap_functions import common_gap_functions
        if name in common_gap_functions:
            f = make_GapElement_Function(self, gap_eval(str(name)))
            assert f.is_function()
            self.__dict__[name] = f
            return f
        else:
            raise AttributeError('No such attribute: '+name+'.')

    def show(self):
        """
        Print statistics about the GAP owned object list

        Slight complication is that we want to do it without accessing
        libgap objects, so we don't create new GapElements as a side
        effect.

        EXAMPLES::

            sage: a = libgap(123)
            sage: b = libgap(456)
            sage: c = libgap(789)
            sage: del b
            sage: libgap.show() # random output
            11 LibGAP elements currently alive
            rec( full := rec( cumulative := 122, deadbags := 9,
            deadkb := 0, freekb := 7785, livebags := 304915,
            livekb := 47367, time := 33, totalkb := 68608 ),
            nfull := 3, npartial := 14 )
        """
        print('{} LibGAP elements currently alive'.format(self.count_GAP_objects()))
        print(self.eval('GasmanStatistics()'))

    def count_GAP_objects(self):
        """
        Return the number of GAP objects that are being tracked by
        libGAP

        OUTPUT:

        An integer

        EXAMPLES::

            sage: libgap.count_GAP_objects()   # random output
            5
        """
        return sum([1 for obj in get_owned_objects()])

    def mem(self):
        """
        Return information about libGAP memory usage

        The GAP workspace is partitioned into 5 pieces (see gasman.c
        in the GAP sources for more details):

        * The **masterpointer area**  contains  all the masterpointers  of  the bags.

        * The **old bags area** contains the bodies of all the  bags that survived at
          least one  garbage collection.  This area is  only  scanned for dead bags
          during a full garbage collection.

        * The **young bags area** contains the bodies of all  the bags that have been
          allocated since the  last garbage collection.  This  area is scanned  for
          dead  bags during  each garbage  collection.

        * The **allocation area** is the storage  that is available for allocation of
          new bags.  When a new bag is allocated the storage for  the body is taken
          from  the beginning of   this area,  and  this  area  is  correspondingly
          reduced.   If  the body does not   fit in the  allocation  area a garbage
          collection is  performed.

        * The **unavailable  area** is  the free  storage that  is not  available for
          allocation.

        OUTPUT:

        This function returns a tuple containing 5 integers. Each is
        the size (in bytes) of the five partitions of the
        workspace. This will potentially change after each GAP garbage
        collection.

        EXAMPLES::

            sage: libgap.collect()
            sage: libgap.mem()   # random output
            (1048576, 6706782, 0, 960930, 0)

            sage: libgap.FreeGroup(3)
            <free group on the generators [ f1, f2, f3 ]>
            sage: libgap.mem()   # random output
            (1048576, 6706782, 47571, 913359, 0)

            sage: libgap.collect()
            sage: libgap.mem()   # random output
            (1048576, 6734785, 0, 998463, 0)
        """
        return memory_usage()

    def collect(self):
        """
        Manually run the garbage collector

        EXAMPLES::

            sage: a = libgap(123)
            sage: del a
            sage: libgap.collect()
        """
        libgap_enter()
        rc = libGAP_CollectBags(0, 1)
        libgap_exit()
        if rc != 1:
            raise RuntimeError('Garbage collection failed.')


libgap = Gap()
