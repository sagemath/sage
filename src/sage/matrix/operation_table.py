r"""
Operation Tables

This module implements general operation tables, which are very matrix-like.
"""

# ****************************************************************************
#       Copyright (C) 2010 Rob Beezer <beezer at ups.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy

from sage.structure.sage_object import SageObject
from sage.matrix.constructor import Matrix


class OperationTable(SageObject):
    r"""
    An object that represents a binary operation as a table.

    Primarily this object is used to provide a
    :meth:`~sage.categories.magmas.Magmas.ParentMethods.multiplication_table`
    for objects in the category of magmas (monoids, groups, ...) and
    :meth:`~sage.categories.additive_magmas.AdditiveMagmas.ParentMethods.addition_table`
    for objects in the category of commutative additive magmas
    (additive monoids, groups, ...).

    INPUT:

    - ``S`` -- a finite algebraic structure (or finite iterable)

    - ``operation`` -- a function of two variables that accepts pairs
        of elements from ``S``. A natural source of such functions is
        the Python :mod:`operator` module, and in particular
        :func:`operator.add` and :func:`operator.mul`. This may also
        be a function defined with ``lambda`` or ``def.``

    - ``names`` -- (default: ``'letters'``)  The type of names
      used, values are:

      * ``'letters'`` -- lowercase ASCII letters are used
        for a base 26 representation of the elements'
        positions in the list given by
        :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
        padded to a common width with leading 'a's.
      * ``'digits'`` -- base 10 representation of the
        elements' positions in the list given by
        :meth:`~sage.matrix.operation_table.OperationTable.column_keys`,
        padded to a common width with leading zeros.
      * ``'elements'`` -- the string representations
        of the elements themselves.
      * a list - a list of strings, where the length
        of the list equals the number of elements.

    - ``elements`` -- (default: ``None``)  A list of elements of ``S``,
      in forms that can be coerced into the structure, eg. their
      string representations. This may be used to impose an alternate
      ordering on the elements of ``S``, perhaps when this is used in
      the context of a particular structure. The default is to use
      whatever ordering the ``S.list()`` method returns. ``elements``
      can also be a subset which is closed under the operation, useful
      perhaps when the set is infinite.

    OUTPUT:

    An object with methods that abstracts multiplication tables,
    addition tables, Cayley tables, etc. It should be general
    enough to be useful for any finite algebraic structure
    whose elements can be combined with a binary operation.
    This is not necessarily meant be constructed directly, but instead
    should be useful for constructing operation tables of various
    algebraic structures that have binary operations.

    EXAMPLES:

    In its most basic use, the table needs a structure and an operation::

        sage: from sage.matrix.operation_table import OperationTable
        sage: G = SymmetricGroup(3)                                                     # needs sage.groups
        sage: OperationTable(G, operation=operator.mul)                                 # needs sage.groups
        *  a b c d e f
         +------------
        a| a b c d e f
        b| b a d c f e
        c| c e a f b d
        d| d f b e a c
        e| e c f a d b
        f| f d e b c a

    With two operations present, we can specify which operation we
    want::

        sage: from sage.matrix.operation_table import OperationTable
        sage: R = Integers(6)
        sage: OperationTable(R, operation=operator.add)
        +  a b c d e f
         +------------
        a| a b c d e f
        b| b c d e f a
        c| c d e f a b
        d| d e f a b c
        e| e f a b c d
        f| f a b c d e

    The default symbol set for elements is lowercase ASCII letters,
    which take on a base 26 flavor for structures with more than
    26 elements. ::

        sage: from sage.matrix.operation_table import OperationTable
        sage: G = DihedralGroup(14)                                                     # needs sage.groups
        sage: OperationTable(G, operator.mul, names='letters')                          # needs sage.groups
         *  aa ab ac ad ae af ag ah ai aj ak al am an ao ap aq ar as at au av aw ax ay az ba bb
          +------------------------------------------------------------------------------------
        aa| aa ab ac ad ae af ag ah ai aj ak al am an ao ap aq ar as at au av aw ax ay az ba bb
        ab| ab aa ad ac af ae ah ag aj ai al ak an am ap ao ar aq at as av au ax aw az ay bb ba
        ac| ac ba aa ae ad ag af ai ah ak aj am al ao an aq ap as ar au at aw av ay ax bb ab az
        ad| ad bb ab af ac ah ae aj ag al ai an ak ap am ar ao at aq av as ax au az aw ba aa ay
        ae| ae az ba ag aa ai ad ak af am ah ao aj aq al as an au ap aw ar ay at bb av ab ac ax
        af| af ay bb ah ab aj ac al ae an ag ap ai ar ak at am av ao ax aq az as ba au aa ad aw
        ag| ag ax az ai ba ak aa am ad ao af aq ah as aj au al aw an ay ap bb ar ab at ac ae av
        ah| ah aw ay aj bb al ab an ac ap ae ar ag at ai av ak ax am az ao ba aq aa as ad af au
        ai| ai av ax ak az am ba ao aa aq ad as af au ah aw aj ay al bb an ab ap ac ar ae ag at
        aj| aj au aw al ay an bb ap ab ar ac at ae av ag ax ai az ak ba am aa ao ad aq af ah as
        ak| ak at av am ax ao az aq ba as aa au ad aw af ay ah bb aj ab al ac an ae ap ag ai ar
        al| al as au an aw ap ay ar bb at ab av ac ax ae az ag ba ai aa ak ad am af ao ah aj aq
        am| am ar at ao av aq ax as az au ba aw aa ay ad bb af ab ah ac aj ae al ag an ai ak ap
        an| an aq as ap au ar aw at ay av bb ax ab az ac ba ae aa ag ad ai af ak ah am aj al ao
        ao| ao ap ar aq at as av au ax aw az ay ba bb aa ab ad ac af ae ah ag aj ai al ak am an
        ap| ap ao aq ar as at au av aw ax ay az bb ba ab aa ac ad ae af ag ah ai aj ak al an am
        aq| aq an ap as ar au at aw av ay ax bb az ab ba ac aa ae ad ag af ai ah ak aj am ao al
        ar| ar am ao at aq av as ax au az aw ba ay aa bb ad ab af ac ah ae aj ag al ai an ap ak
        as| as al an au ap aw ar ay at bb av ab ax ac az ae ba ag aa ai ad ak af am ah ao aq aj
        at| at ak am av ao ax aq az as ba au aa aw ad ay af bb ah ab aj ac al ae an ag ap ar ai
        au| au aj al aw an ay ap bb ar ab at ac av ae ax ag az ai ba ak aa am ad ao af aq as ah
        av| av ai ak ax am az ao ba aq aa as ad au af aw ah ay aj bb al ab an ac ap ae ar at ag
        aw| aw ah aj ay al bb an ab ap ac ar ae at ag av ai ax ak az am ba ao aa aq ad as au af
        ax| ax ag ai az ak ba am aa ao ad aq af as ah au aj aw al ay an bb ap ab ar ac at av ae
        ay| ay af ah bb aj ab al ac an ae ap ag ar ai at ak av am ax ao az aq ba as aa au aw ad
        az| az ae ag ba ai aa ak ad am af ao ah aq aj as al au an aw ap ay ar bb at ab av ax ac
        ba| ba ac ae aa ag ad ai af ak ah am aj ao al aq an as ap au ar aw at ay av bb ax az ab
        bb| bb ad af ab ah ac aj ae al ag an ai ap ak ar am at ao av aq ax as az au ba aw ay aa

    Another symbol set is base 10 digits, padded with leading
    zeros to make a common width. ::

        sage: from sage.matrix.operation_table import OperationTable
        sage: G = AlternatingGroup(4)                                                   # needs sage.groups
        sage: OperationTable(G, operator.mul, names='digits')                           # needs sage.groups
         *  00 01 02 03 04 05 06 07 08 09 10 11
          +------------------------------------
        00| 00 01 02 03 04 05 06 07 08 09 10 11
        01| 01 02 00 05 03 04 07 08 06 11 09 10
        02| 02 00 01 04 05 03 08 06 07 10 11 09
        03| 03 06 09 00 07 10 01 04 11 02 05 08
        04| 04 08 10 02 06 11 00 05 09 01 03 07
        05| 05 07 11 01 08 09 02 03 10 00 04 06
        06| 06 09 03 10 00 07 04 11 01 08 02 05
        07| 07 11 05 09 01 08 03 10 02 06 00 04
        08| 08 10 04 11 02 06 05 09 00 07 01 03
        09| 09 03 06 07 10 00 11 01 04 05 08 02
        10| 10 04 08 06 11 02 09 00 05 03 07 01
        11| 11 05 07 08 09 01 10 02 03 04 06 00

    If the group's elements are not too cumbersome,
    or the group is small, then the string representation
    of the elements can be used. ::

        sage: from sage.matrix.operation_table import OperationTable
        sage: G = AlternatingGroup(3)                                                   # needs sage.groups
        sage: OperationTable(G, operator.mul, names='elements')                         # needs sage.groups
              *       () (1,2,3) (1,3,2)
               +------------------------
             ()|      () (1,2,3) (1,3,2)
        (1,2,3)| (1,2,3) (1,3,2)      ()
        (1,3,2)| (1,3,2)      () (1,2,3)

    You can give the elements any names you like, but they need to be ordered
    in the same order as returned by the
    :meth:`~sage.matrix.operation_table.OperationTable.column_keys`
    method.  ::

        sage: # needs sage.groups
        sage: from sage.matrix.operation_table import OperationTable
        sage: G = QuaternionGroup()
        sage: T = OperationTable(G, operator.mul)
        sage: T.column_keys()
        ((), (1,2,3,4)(5,6,7,8), ..., (1,8,3,6)(2,7,4,5))
        sage: names=['1', 'I', '-1', '-I', 'J', '-K', '-J', 'K']
        sage: T.change_names(names=names)
        sage: sorted(T.translation().items())
        [('-1', (1,3)(2,4)(5,7)(6,8)), ..., ('K', (1,8,3,6)(2,7,4,5))]
        sage: T
         *   1  I -1 -I  J -K -J  K
          +------------------------
         1|  1  I -1 -I  J -K -J  K
         I|  I -1 -I  1  K  J -K -J
        -1| -1 -I  1  I -J  K  J -K
        -I| -I  1  I -1 -K -J  K  J
         J|  J -K -J  K -1 -I  1  I
        -K| -K -J  K  J  I -1 -I  1
        -J| -J  K  J -K  1  I -1 -I
         K|  K  J -K -J -I  1  I -1

    With the right functions and a list comprehension, custom
    names can be easier.  A multiplication table for hex digits
    (without carries)::

        sage: from sage.matrix.operation_table import OperationTable
        sage: R = Integers(16)
        sage: names=['{:x}'.format(Integer(a)) for a in R]
        sage: OperationTable(R, operation=operator.mul, names=names)
        *  0 1 2 3 4 5 6 7 8 9 a b c d e f
         +--------------------------------
        0| 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        1| 0 1 2 3 4 5 6 7 8 9 a b c d e f
        2| 0 2 4 6 8 a c e 0 2 4 6 8 a c e
        3| 0 3 6 9 c f 2 5 8 b e 1 4 7 a d
        4| 0 4 8 c 0 4 8 c 0 4 8 c 0 4 8 c
        5| 0 5 a f 4 9 e 3 8 d 2 7 c 1 6 b
        6| 0 6 c 2 8 e 4 a 0 6 c 2 8 e 4 a
        7| 0 7 e 5 c 3 a 1 8 f 6 d 4 b 2 9
        8| 0 8 0 8 0 8 0 8 0 8 0 8 0 8 0 8
        9| 0 9 2 b 4 d 6 f 8 1 a 3 c 5 e 7
        a| 0 a 4 e 8 2 c 6 0 a 4 e 8 2 c 6
        b| 0 b 6 1 c 7 2 d 8 3 e 9 4 f a 5
        c| 0 c 8 4 0 c 8 4 0 c 8 4 0 c 8 4
        d| 0 d a 7 4 1 e b 8 5 2 f c 9 6 3
        e| 0 e c a 8 6 4 2 0 e c a 8 6 4 2
        f| 0 f e d c b a 9 8 7 6 5 4 3 2 1

    This should be flexible enough to create a variety
    of such tables.  ::

        sage: from sage.matrix.operation_table import OperationTable
        sage: from operator import xor
        sage: T=OperationTable(ZZ, xor, elements=range(8))
        sage: T
        .  a b c d e f g h
         +----------------
        a| a b c d e f g h
        b| b a d c f e h g
        c| c d a b g h e f
        d| d c b a h g f e
        e| e f g h a b c d
        f| f e h g b a d c
        g| g h e f c d a b
        h| h g f e d c b a
        sage: names=['000', '001','010','011','100','101','110','111']
        sage: T.change_names(names)
        sage: T.set_print_symbols('^', '\\land')
        sage: T
          ^  000 001 010 011 100 101 110 111
           +--------------------------------
        000| 000 001 010 011 100 101 110 111
        001| 001 000 011 010 101 100 111 110
        010| 010 011 000 001 110 111 100 101
        011| 011 010 001 000 111 110 101 100
        100| 100 101 110 111 000 001 010 011
        101| 101 100 111 110 001 000 011 010
        110| 110 111 100 101 010 011 000 001
        111| 111 110 101 100 011 010 001 000

        sage: T = OperationTable([False, True], operator.or_, names = 'elements')
        sage: T
            .  False  True
             +------------
        False| False  True
         True|  True  True


    TESTS:

    Empty structures behave acceptably, though the ASCII table looks a bit
    odd.  The LaTeX version works much better. ::

        sage: from sage.matrix.operation_table import OperationTable
        sage: L = FiniteSemigroups().example(())
        sage: L
        An example of a finite semigroup: the left regular band generated by ()
        sage: T = OperationTable(L, operation=operator.mul)
        sage: T
        *
         +
        sage: T._latex_()
        '{\\setlength{\\arraycolsep}{2ex}\n\\begin{array}{r|*{0}{r}}\n\\multicolumn{1}{c|}{\\ast}\\\\\\hline\n\\end{array}}'

    If the algebraic structure cannot be listed (like when it is infinite)
    then there is no way to create a table. ::

        sage: from sage.matrix.operation_table import OperationTable
        sage: OperationTable(ZZ, operator.mul)
        Traceback (most recent call last):
        ...
        ValueError: Integer Ring is infinite

    The value of ``elements`` must be a subset of the algebraic
    structure, in forms that can be coerced into the structure.
    Here we demonstrate the proper use first::

        sage: # needs sage.groups
        sage: from sage.matrix.operation_table import OperationTable
        sage: H = CyclicPermutationGroup(4)
        sage: H.list()
        [(), (1,2,3,4), (1,3)(2,4), (1,4,3,2)]
        sage: elts = ['()', '(1,3)(2,4)']
        sage: OperationTable(H, operator.mul, elements=elts)
        *  a b
         +----
        a| a b
        b| b a

    This can be rewritten so as to pass the actual elements of the
    group ``H``, using a simple ``for`` loop::

        sage: # needs sage.groups
        sage: L = H.list()    #list of elements of the group H
        sage: elts = [L[i] for i in {0, 2}]
        sage: elts
        [(), (1,3)(2,4)]
        sage: OperationTable(H, operator.mul, elements=elts)
        *  a b
         +----
        a| a b
        b| b a

    Here are a couple of improper uses::

        sage: # needs sage.groups
        sage: elts.append(5)
        sage: OperationTable(H, operator.mul, elements=elts)
        Traceback (most recent call last):
        ...
        TypeError: unable to coerce 5 into Cyclic group of order 4 as a permutation group
        sage: elts[2] = '(1,3,2,4)'
        sage: OperationTable(H, operator.mul, elements=elts)
        Traceback (most recent call last):
        ...
        TypeError: unable to coerce (1,3,2,4) into Cyclic group of order 4 as a permutation group
        sage: elts[2] = '(1,2,3,4)'
        sage: OperationTable(H, operator.mul, elements=elts)
        Traceback (most recent call last):
        ...
        ValueError: (1,3)(2,4)*(1,2,3,4)=(1,4,3,2), and so the set is not closed

    Unusable functions should be recognized as such::

        sage: H = CyclicPermutationGroup(4)                                             # needs sage.groups
        sage: OperationTable(H, operator.add)                                           # needs sage.groups
        Traceback (most recent call last):
        ...
        TypeError: elements () and () of Cyclic group of order 4 as a permutation group are incompatible with operation: <built-in function add>
        sage: from operator import xor
        sage: OperationTable(H, xor)                                                    # needs sage.groups
        Traceback (most recent call last):
        ...
        TypeError: elements () and () of Cyclic group of order 4 as a permutation group are incompatible with operation: <built-in function xor>

    We construct the multiplication table for a finite finitely presented
    group, where there is no normalization done when computing the hash::

        sage: GU.<s,t> = FreeGroup()                                                    # needs sage.groups
        sage: gr0 = GU / (s^(-2)*t*s*t, t^(-2)*s*t*s, s*t*s*t)                          # needs sage.groups
        sage: gr0.multiplication_table()                                                # needs sage.groups
        *  a b c d e f g h i j k l
         +------------------------
        a| a b c d e f g h i j k l
        b| b e f g a i j k c d l h
        c| c g a h l k b d j i f e
        d| d k h a i g f c e l b j
        e| e a i j b c d l f g h k
        f| f j b k h l e g d c i a
        g| g l k b c j i f a h e d
        h| h f d c j b k a l e g i
        i| i d e l k h a j g f c b
        j| j h l e f d c i b k a g
        k| k i g f d e l b h a j c
        l| l c j i g a h e k b d f

    AUTHORS:

    - Rob Beezer (2010-03-15)
    - Bruno Edwards (2022-10-31)
    """

    def __init__(self, S, operation, names='letters', elements=None):
        r"""
        TESTS::

            sage: from sage.matrix.operation_table import OperationTable
            sage: G = SymmetricGroup(3)                                                 # needs sage.groups
            sage: T = OperationTable(G, operator.mul)                                   # needs sage.groups
            sage: TestSuite(T).run()                                                    # needs sage.groups
        """
        # Determine the elements of S, specified or not
        # If elements are given, we check if they are all in S
        # Note: there exist listable infinite objects (like ZZ)
        if (elements is None):
            if hasattr(S, 'is_finite'):
                if not S.is_finite():
                    raise ValueError('%s is infinite' % S)
            try:
                try:
                    elems = sorted(S)
                except TypeError:
                    elems = list(S)
            except Exception:
                raise ValueError('unable to determine elements of %s' % S)
        else:
            elems = []
            try:
                for e in elements:
                    coerced = S(e)
                    if coerced not in elems:
                        elems.append(coerced)
            except Exception:
                raise TypeError('unable to coerce %s into %s' % (e, S))
        self._elts = tuple(elems)
        self._n = len(self._elts)
        self._name_dict = {}

        # Map elements to strings
        self._width, self._names, self._name_dict = self._name_maker(names)

        # Determine the operation, if given by a string
        # Some simple symbols are supported,
        # Add support for other symbols in the supported dictionary
        # Triple is (function, ascii-symbol, latex-symbol)
        # ascii-symbol must be exactly one character wide
        # Note double-backslash to escape properly for latex
        from operator import add, mul
        supported = {
            add: (add, '+', '+'),
            mul: (mul, '*', '\\ast')
        }
        # default symbols for upper-left-hand-corner of table
        self._ascii_symbol = '.'
        self._latex_symbol = '\\cdot'
        if operation in supported.keys():
            chosen = supported[operation]
            operation = chosen[0]
            self._ascii_symbol = chosen[1]
            self._latex_symbol = chosen[2]
        self._operation = operation
        # We assume now that operation is a function that acts
        # as a closed binary operation on the elements.
        # If not, we'll discover that next in actual use.

        self._table = []

        # the elements might not be hashable. But if they are it is much
        # faster to lookup in a hash table rather than in a list!
        try:
            get_row = {e: i for i, e in enumerate(self._elts)}.__getitem__
        except TypeError:
            get_row = self._elts.index

        for g in self._elts:
            row = []
            for h in self._elts:
                try:
                    result = self._operation(g, h)
                except Exception:
                    raise TypeError('elements %s and %s of %s are incompatible with operation: %s' % (
                        g, h, S, self._operation))

                try:
                    r = get_row(result)
                except (KeyError, ValueError):
                    failed = True
                    # There might be an issue with the hashing, fall back to
                    #   getting the index (which simply uses ==).
                    if get_row != self._elts.index:
                        failed = False
                        get_row = self._elts.index
                        try:
                            r = get_row(result)
                        except (KeyError, ValueError):
                            failed = True
                    if failed:
                        raise ValueError('%s%s%s=%s, and so the set is not closed' % (
                            g, self._ascii_symbol, h, result))

                row.append(r)
            self._table.append(row)

    def _name_maker(self, names):
        r"""
        Helper function to create names of elements of algebraic structures.

        INPUT:

        Identical to the input for :class:`OperationTable` and :meth:`change_names`,
        so look there for details.

        OUTPUT:

        - ``width`` -- integer giving the maximum width of the strings
          describing the elements.  This is used for formatting the ASCII
          version of the table.
        - ``name_list`` -- list of strings naming the elements, in the
          same order as given by the :meth:`list` method
        - ``name_dict`` -- dictionary giving the correspondence between the
          strings and the actual elements.  So the keys are the strings and
          the values are the elements of the structure.

        EXAMPLES:

        This routine is tested extensively in the :class:`OperationTable`
        and :meth:`change_names` methods.  So we just demonstrate
        the nature of the output here. ::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = SymmetricGroup(3)
            sage: T = OperationTable(G, operator.mul)
            sage: w, l, d = T._name_maker('letters')
            sage: w
            1
            sage: l[0]
            'a'
            sage: d['a']
            ()

        TESTS:

        We test the error conditions here, rather than as part of the
        doctests for the :class:`OperationTable` and :meth:`change_names`
        methods that rely on this one. ::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = AlternatingGroup(3)
            sage: T = OperationTable(G, operator.mul)
            sage: T._name_maker(['x'])
            Traceback (most recent call last):
            ...
            ValueError: list of element names must be the same size as the set, 1 != 3
            sage: T._name_maker(['x', 'y', 4])
            Traceback (most recent call last):
            ...
            ValueError: list of element names must only contain strings, not 4
            sage: T._name_maker('blatzo')
            Traceback (most recent call last):
            ...
            ValueError: element names must be a list, or one of the keywords: 'letters', 'digits', 'elements'
        """
        from math import log, log10
        name_list = []
        if names == 'digits':
            if self._n == 0 or self._n == 1:
                width = 1
            else:
                width = int(log10(self._n - 1)) + 1
            for i in range(self._n):
                name_list.append('{0:0{1}d}'.format(i, width))
        elif names == 'letters':
            from string import ascii_lowercase as letters
            from sage.rings.integer import Integer
            base = len(letters)
            if self._n == 0 or self._n == 1:
                width = 1
            else:
                width = int(log(self._n - 1, base)) + 1
            for i in range(self._n):
                places = Integer(i).digits(
                    base=base, digits=letters, padto=width)
                places.reverse()
                name_list.append(''.join(places))
        elif names == 'elements':
            width = 0
            for e in self._elts:
                estr = repr(e)
                width = max(len(estr), width)
                name_list.append(estr)
        elif isinstance(names, list):
            if len(names) != self._n:
                raise ValueError('list of element names must be the same size as the set, %s != %s' % (
                    len(names), self._n))
            width = 0
            for name in names:
                if not isinstance(name, str):
                    raise ValueError(
                        'list of element names must only contain strings, not %s' % name)
                width = max(len(name), width)
                name_list.append(name)
        else:
            raise ValueError(
                "element names must be a list, or one of the keywords: 'letters', 'digits', 'elements'")
        name_dict = {}
        for i in range(self._n):
            name_dict[name_list[i]] = self._elts[i]
        return width, name_list, name_dict

    def __getitem__(self, pair):
        r"""
        Return the element of the table, given the elements indexing its position.

        INPUT:

        - ``pair`` -- two elements of the structure

        OUTPUT:

        The element of the structure computed by the operation for
        the two input elements (in the order provided).

        This uses the table as a look-up device.  If you want to use
        the operation, then use the operation.

        EXAMPLES::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = DiCyclicGroup(3)
            sage: T = OperationTable(G, operator.mul)
            sage: T.column_keys()
            ((), (5,6,7), ..., (1,4,2,3)(5,7))
            sage: T[G('(1,2)(3,4)(5,6,7)'), G('(1,3,2,4)(5,7)')]
            (1,4,2,3)(5,6)

        TESTS::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = DiCyclicGroup(3)
            sage: T = OperationTable(G, operator.mul)
            sage: T[G('(1,2)(3,4)(5,6,7)')]
            Traceback (most recent call last):
            ...
            TypeError: indexing into an operation table requires exactly two elements
            sage: T[G('(1,2)(3,4)(5,6,7)'), G('(1,3,2,4)(5,7)'), G('(1,3,2,4)(5,7)')]
            Traceback (most recent call last):
            ...
            TypeError: indexing into an operation table requires exactly two elements
            sage: T[2, 3]
            Traceback (most recent call last):
            ...
            IndexError: invalid indices of operation table: (2, 3)
            sage: T['(1,512)', '(1,3,2,4)(5,7)']
            Traceback (most recent call last):
            ...
            IndexError: invalid indices of operation table: ((1,512), (1,3,2,4)(5,7))
        """
        if not (isinstance(pair, tuple) and len(pair) == 2):
            raise TypeError(
                'indexing into an operation table requires exactly two elements')
        g, h = pair
        try:
            row = self._elts.index(g)
            col = self._elts.index(h)
        except ValueError:
            raise IndexError(
                'invalid indices of operation table: (%s, %s)' % (g, h))
        return self._elts[self._table[row][col]]

    def __eq__(self, other):
        r"""
        Return the comparison between two tables.

        INPUT:

        - ``other`` -- a second table to compare to ``self``

        OUTPUT:

        Tables are equal if they have the same operation and elements.

        EXAMPLES::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = CyclicPermutationGroup(6)
            sage: H = CyclicPermutationGroup(3)
            sage: P = OperationTable(G, operator.mul)
            sage: Q = OperationTable(G, operator.mul)
            sage: R = OperationTable(H, operator.mul)
            sage: S = OperationTable(G, operator.truediv)
            sage: P == P, P == Q, P == R, P == S
            (True, True, False, False)
        """
        return (self._elts == other._elts) and (self._operation == other._operation)

    def __ne__(self, other):
        """
        Inequality test, by negation of :meth:`.__eq__`.

        EXAMPLES::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = CyclicPermutationGroup(6)
            sage: H = CyclicPermutationGroup(3)
            sage: P = OperationTable(G, operator.mul)
            sage: Q = OperationTable(G, operator.mul)
            sage: R = OperationTable(H, operator.mul)
            sage: S = OperationTable(G, operator.truediv)
            sage: P != P, P != Q, P != R, P != S
            (False, False, True, True)
        """
        return not self == other

    def _repr_(self):
        r"""
        Return a printable version of the operation table.

        EXAMPLES::

            sage: from sage.matrix.operation_table import OperationTable
            sage: R = Integers(5)
            sage: T = OperationTable(R, operation=operator.add)
            sage: print(T._repr_())
            +  a b c d e
             +----------
            a| a b c d e
            b| b c d e a
            c| c d e a b
            d| d e a b c
            e| e a b c d
        """
        return self._ascii_table()

    def set_print_symbols(self, ascii, latex):
        r"""
        Set the symbols used for text and LaTeX printing of operation tables.

        INPUT:

        - ``ascii`` -- a single character for text table
        - ``latex`` -- string to represent an operation in LaTeX math mode;
          note the need for double-backslashes to escape properly

        EXAMPLES::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = AlternatingGroup(3)
            sage: T = OperationTable(G, operator.mul)
            sage: T.set_print_symbols('@', '\\times')
            sage: T
            @  a b c
             +------
            a| a b c
            b| b c a
            c| c a b
            sage: T._latex_()
            '{\\setlength{\\arraycolsep}{2ex}\n\\begin{array}{r|*{3}{r}}\n\\multicolumn{1}{c|}{\\times}&a&b&c\\\\\\hline\n{}a&a&b&c\\\\\n{}b&b&c&a\\\\\n{}c&c&a&b\\\\\n\\end{array}}'

        TESTS::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = AlternatingGroup(3)
            sage: T = OperationTable(G, operator.mul)
            sage: T.set_print_symbols('@', 5)
            Traceback (most recent call last):
            ...
            ValueError: LaTeX symbol must be a string, not 5
            sage: T.set_print_symbols('@x@', '\\times')
            Traceback (most recent call last):
            ...
            ValueError: ASCII symbol should be a single character, not @x@
            sage: T.set_print_symbols(5, '\\times')
            Traceback (most recent call last):
            ...
            ValueError: ASCII symbol should be a single character, not 5
        """
        if not isinstance(ascii, str) or not len(ascii) == 1:
            raise ValueError(
                'ASCII symbol should be a single character, not %s' % ascii)
        if not isinstance(latex, str):
            raise ValueError('LaTeX symbol must be a string, not %s' % latex)
        self._ascii_symbol = ascii
        self._latex_symbol = latex
        return None

    def column_keys(self):
        r"""
        Return a tuple of the elements used to build the table.

        .. NOTE:: ``column_keys`` and ``row_keys`` are identical.
           Both list the elements in the order used to label the table.

        OUTPUT:

        The elements of the algebraic structure used to build
        the table, as a list.  But most importantly, elements are
        present in the list in the order which they appear in
        the table's column headings.

        EXAMPLES::

            sage: from sage.matrix.operation_table import OperationTable
            sage: G = AlternatingGroup(3)                                               # needs sage.groups
            sage: T = OperationTable(G, operator.mul)                                   # needs sage.groups
            sage: T.column_keys()                                                       # needs sage.groups
            ((), (1,2,3), (1,3,2))
        """
        return self._elts

    # The ordered list of row and column elements are identical
    # given the current design, so these methods are aliases.  If
    # expanded to allow different orderings (maybe interesting in
    # non-commutative cases?), then these will need to be
    # implemented separately.
    row_keys = column_keys

    def translation(self):
        r"""
        Return a dictionary associating names with elements.

        OUTPUT:

        A dictionary whose keys are strings used as names
        for entries of the table and values that are the
        actual elements of the algebraic structure.

        EXAMPLES::

            sage: from sage.matrix.operation_table import OperationTable
            sage: G = AlternatingGroup(3)                                               # needs sage.groups
            sage: T = OperationTable(G, operator.mul, names=['p','q','r'])              # needs sage.groups
            sage: T.translation()                                                       # needs sage.groups
            {'p': (), 'q': (1,2,3), 'r': (1,3,2)}
        """
        return self._name_dict

    def table(self):
        r"""
        Return the table as a list of lists,
        using integers to reference the elements.

        OUTPUT:

        The rows of the table, as a list of rows, each row
        being a list of integer entries.  The integers correspond
        to the order of the elements in the headings of the table
        and the order of the output of the :meth:`list` method.

        EXAMPLES::

            sage: from sage.matrix.operation_table import OperationTable
            sage: C = CyclicPermutationGroup(3)                                         # needs sage.groups
            sage: T=OperationTable(C, operator.mul)                                     # needs sage.groups
            sage: T.table()                                                             # needs sage.groups
            [[0, 1, 2], [1, 2, 0], [2, 0, 1]]
        """
        return self._table

    def change_names(self, names):
        r"""
        For an existing operation table, change the names used for the elements.

        INPUT:

        - ``names`` -- the type of names used, values are:

          * ``'letters'`` -- lowercase ASCII letters are used
            for a base 26 representation of the elements'
            positions in the list given by :meth:`list`,
            padded to a common width with leading 'a's.
          * ``'digits'`` -- base 10 representation of the
            elements' positions in the list given by
            :meth:`list`, padded to a common width
            with leading zeros.
          * ``'elements'`` -- the string representations
            of the elements themselves.
          * a list - a list of strings, where the length
            of the list equals the number of elements.

        OUTPUT:
        ``None``.  This method changes the table "in-place",
        so any printed version will change and the output of
        the :meth:`dict` will also change.  So any items of
        interest about a particular table need to be copied/saved
        prior to calling this method.

        EXAMPLES:

        More examples can be found in the documentation for
        :class:`OperationTable` since creating a new
        operation table uses the same routine. ::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: D = DihedralGroup(2)
            sage: T = OperationTable(D, operator.mul)
            sage: T
            *  a b c d
             +--------
            a| a b c d
            b| b a d c
            c| c d a b
            d| d c b a
            sage: T.translation()['c']
            (1,2)
            sage: T.change_names('digits')
            sage: T
            *  0 1 2 3
             +--------
            0| 0 1 2 3
            1| 1 0 3 2
            2| 2 3 0 1
            3| 3 2 1 0
            sage: T.translation()['2']
            (1,2)
            sage: T.change_names('elements')
            sage: T
                     *          ()      (3,4)      (1,2) (1,2)(3,4)
                      +--------------------------------------------
                    ()|         ()      (3,4)      (1,2) (1,2)(3,4)
                 (3,4)|      (3,4)         () (1,2)(3,4)      (1,2)
                 (1,2)|      (1,2) (1,2)(3,4)         ()      (3,4)
            (1,2)(3,4)| (1,2)(3,4)      (1,2)      (3,4)         ()
            sage: T.translation()['(1,2)']
            (1,2)
            sage: T.change_names(['w', 'x', 'y', 'z'])
            sage: T
            *  w x y z
             +--------
            w| w x y z
            x| x w z y
            y| y z w x
            z| z y x w
            sage: T.translation()['y']
            (1,2)
        """
        self._width, self._names, self._name_dict = self._name_maker(names)
        return None

    def matrix_of_variables(self):
        r"""
        This method provides some backward compatibility for
        Cayley tables of groups, whose output was
        restricted to this single format.

        EXAMPLES:

        The output here is from the doctests for the old
        ``cayley_table()`` method for permutation groups. ::

            sage: # needs sage.groups
            sage: from sage.matrix.operation_table import OperationTable
            sage: G = PermutationGroup(['(1,2,3)', '(2,3)'])
            sage: T = OperationTable(G, operator.mul)
            sage: T.matrix_of_variables()
            [x0 x1 x2 x3 x4 x5]
            [x1 x0 x3 x2 x5 x4]
            [x2 x4 x0 x5 x1 x3]
            [x3 x5 x1 x4 x0 x2]
            [x4 x2 x5 x0 x3 x1]
            [x5 x3 x4 x1 x2 x0]
            sage: T.column_keys()[2]*T.column_keys()[2] == T.column_keys()[0]
            True
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.matrix.matrix_space import MatrixSpace
        from sage.rings.rational_field import QQ
        R = PolynomialRing(QQ, 'x', self._n)
        MS = MatrixSpace(R, self._n, self._n)
        entries = [R('x'+str(self._table[i][j]))
                   for i in range(self._n) for j in range(self._n)]
        return MS(entries)

    def color_table(self, element_names=True, cmap=None, **options):
        r"""
        Return a graphic image as a square grid where entries are color coded.

        INPUT:

        - ``element_names`` -- (default: ``True``) whether to display text with
          element names on the image

        - ``cmap`` -- (default: :obj:`matplotlib.cm.gist_rainbow`) color map for plot, see :mod:`matplotlib.cm`

        - ``**options`` -- passed on to :func:`~sage.plot.matrix_plot.matrix_plot`

        EXAMPLES::

            sage: from sage.matrix.operation_table import OperationTable
            sage: OTa = OperationTable(SymmetricGroup(3), operation=operator.mul)       # needs sage.groups
            sage: OTa.color_table()                                                     # needs sage.groups sage.plot
            Graphics object consisting of 37 graphics primitives

        .. PLOT::

            from sage.matrix.operation_table import OperationTable
            OTa = OperationTable(SymmetricGroup(3), operation=operator.mul)
            sphinx_plot(OTa.color_table(), figsize=(3.0, 3.0))
        """
        from sage.plot.matrix_plot import matrix_plot
        from sage.plot.text import text

        if cmap is None:
            from matplotlib.cm import gist_rainbow as cmap

        # Base matrix plot object, without text
        plot = matrix_plot(Matrix(self._table), cmap=cmap,
                           frame=False, **options)

        if element_names:

            # adapted from ._ascii_table()
            # prepare widenames[] list for labelling on image
            n = self._n
            width = self._width

            widenames = []
            for name in self._names:
                widenames.append("{0: >{1}s}".format(name, width))

            # iterate through each element
            for g in range(n):
                for h in range(n):
                    # add text to the plot
                    tPos = (h, g)
                    tText = widenames[self._table[g][h]]
                    t = text(tText, tPos, rgbcolor=(0, 0, 0))
                    plot = plot + t

        # https://moyix.blogspot.com/2022/09/someones-been-messing-with-my-subnormals.html
        import warnings
        warnings.filterwarnings("ignore", message="The value of the smallest subnormal for")

        return plot

    def gray_table(self, **options):
        r"""
        Return a graphic image as a square grid where entries are displayed in grayscale.

        INPUT:

        - ``element_names`` -- boolean (default: ``True``); whether to display
          text with element names on the image

        - ``**options`` -- passed on to :func:`~sage.plot.matrix_plot.matrix_plot`

        EXAMPLES::

            sage: from sage.matrix.operation_table import OperationTable
            sage: OTa = OperationTable(SymmetricGroup(3), operation=operator.mul)       # needs sage.groups
            sage: OTa.gray_table()                                                      # needs sage.groups sage.plot
            Graphics object consisting of 37 graphics primitives

        .. PLOT::

            from sage.matrix.operation_table import OperationTable
            OTa = OperationTable(SymmetricGroup(3), operation=operator.mul)
            sphinx_plot(OTa.gray_table(), figsize=(3.0, 3.0))
        """
        from matplotlib.cm import Greys
        return self.color_table(cmap=Greys, **options)

    def _ascii_table(self):
        r"""
        Return a string that is an ASCII version of the table.

        EXAMPLES::

            sage: from sage.matrix.operation_table import OperationTable
            sage: R = Integers(5)
            sage: T = OperationTable(R, operator.add)
            sage: print(T._ascii_table())
            +  a b c d e
             +----------
            a| a b c d e
            b| b c d e a
            c| c d e a b
            d| d e a b c
            e| e a b c d

        The table should adjust its column width to accommodate the width of the
        strings used to represent elements.  ::

            sage: from sage.matrix.operation_table import OperationTable
            sage: R = Integers(10)
            sage: T = OperationTable(R, operator.mul, names='digits')
            sage: print(T._ascii_table())
            *  0 1 2 3 4 5 6 7 8 9
             +--------------------
            0| 0 0 0 0 0 0 0 0 0 0
            1| 0 1 2 3 4 5 6 7 8 9
            2| 0 2 4 6 8 0 2 4 6 8
            3| 0 3 6 9 2 5 8 1 4 7
            4| 0 4 8 2 6 0 4 8 2 6
            5| 0 5 0 5 0 5 0 5 0 5
            6| 0 6 2 8 4 0 6 2 8 4
            7| 0 7 4 1 8 5 2 9 6 3
            8| 0 8 6 4 2 0 8 6 4 2
            9| 0 9 8 7 6 5 4 3 2 1

        ::

            sage: from sage.matrix.operation_table import OperationTable
            sage: R = Integers(11)
            sage: T = OperationTable(R, operator.mul, names='digits')
            sage: print(T._ascii_table())
             *  00 01 02 03 04 05 06 07 08 09 10
              +---------------------------------
            00| 00 00 00 00 00 00 00 00 00 00 00
            01| 00 01 02 03 04 05 06 07 08 09 10
            02| 00 02 04 06 08 10 01 03 05 07 09
            03| 00 03 06 09 01 04 07 10 02 05 08
            04| 00 04 08 01 05 09 02 06 10 03 07
            05| 00 05 10 04 09 03 08 02 07 01 06
            06| 00 06 01 07 02 08 03 09 04 10 05
            07| 00 07 03 10 06 02 09 05 01 08 04
            08| 00 08 05 02 10 07 04 01 09 06 03
            09| 00 09 07 05 03 01 10 08 06 04 02
            10| 00 10 09 08 07 06 05 04 03 02 01

        ::

            sage: from sage.matrix.operation_table import OperationTable
            sage: R = Integers(4)
            sage: T = OperationTable(R, operator.mul, names=['x','y','wwww', 'z'])
            sage: print(T._ascii_table())
               *     x    y wwww    z
                +--------------------
               x|    x    x    x    x
               y|    x    y wwww    z
            wwww|    x wwww    x wwww
               z|    x    z wwww    y
        """
        n = self._n
        width = self._width

        widenames = []
        for name in self._names:
            widenames.append('{0: >{1}s}'.format(name, width))

        # Headers
        table = ['{0: >{1}s} '.format(self._ascii_symbol, width)]
        table += [' '+widenames[i] for i in range(n)]+['\n']
        table += [' ']*width + ['+'] + ['-']*(n*(width+1))+['\n']

        # Row labels, body of table
        for g in range(n):
            table.append(widenames[g]+'|')
            for h in range(n):
                table.append(' '+widenames[self._table[g][h]])
            table.append('\n')
        return ''.join(table)

    def _latex_(self):
        r"""
        Return a `LaTeX` version of the operation table as a string,
        using a `LaTeX` ``array`` environment.

        EXAMPLES::

            sage: from sage.matrix.operation_table import OperationTable
            sage: R = Integers(2)
            sage: T = OperationTable(R, operation=operator.mul)
            sage: T._latex_()
            '{\\setlength{\\arraycolsep}{2ex}\n\\begin{array}{r|*{2}{r}}\n\\multicolumn{1}{c|}{\\ast}&a&b\\\\\\hline\n{}a&a&a\\\\\n{}b&a&b\\\\\n\\end{array}}'
        """
        n = self._n
        names = self._names

        # Headers
        table = ['{\\setlength{\\arraycolsep}{2ex}\n']
        table.append('\\begin{array}{r|*{'+str(n)+'}{r}}\n')
        table.append('\\multicolumn{1}{c|}{'+self._latex_symbol+'}')
        table += ['&'+names[i] for i in range(n)]
        table.append('\\\\\\hline\n')

        # Row label and body of table
        for g in range(n):
            # Interrupts newline and [], so not line spacing
            table.append('{}')
            table.append(names[g])
            for h in range(n):
                table.append('&'+names[self._table[g][h]])
            table.append('\\\\\n')

        # Finish
        table.append('\\end{array}')
        table.append('}')
        return ''.join(table)
