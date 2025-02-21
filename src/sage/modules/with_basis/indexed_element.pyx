r"""
An element in an indexed free module

AUTHORS:

- Travis Scrimshaw (03-2017): Moved code from :mod:`sage.combinat.free_module`.
- Travis Scrimshaw (29-08-2022): Implemented ``copy`` as the identity map.
"""

# ****************************************************************************
#       Copyright (C) 2017, 2022 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element cimport parent
from sage.structure.richcmp cimport richcmp, rich_to_bool
from cpython.object cimport Py_NE, Py_EQ

from sage.misc.repr import repr_lincomb
from sage.misc.superseded import deprecation
from sage.typeset.ascii_art import AsciiArt, empty_ascii_art, ascii_art
from sage.typeset.unicode_art import UnicodeArt, empty_unicode_art, unicode_art
from sage.data_structures.blas_dict cimport add, negate, scal, axpy
from sage.categories.modules import _Fields


cdef class IndexedFreeModuleElement(ModuleElement):
    r"""
    Element class for :class:`~sage.combinat.free_module.CombinatorialFreeModule`.

    TESTS::

        sage: import collections.abc
        sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
        sage: B = F.basis()
        sage: f = B['a'] + 3*B['c']; f
        B['a'] + 3*B['c']
        sage: isinstance(f, collections.abc.Sized)
        True
        sage: isinstance(f, collections.abc.Iterable)
        True
        sage: isinstance(f, collections.abc.Collection)  # known bug - will be fixed by removing __contains__
        False
    """
    def __init__(self, M, x):
        """
        Create a combinatorial module element.

        This should never be called directly, but only through the
        parent combinatorial free module's :meth:`__call__` method.

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']; f
            B['a'] + 3*B['c']
            sage: f == loads(dumps(f))
            True
        """
        ModuleElement.__init__(self, M)
        self._monomial_coefficients = x
        self._hash_set = False

    def __iter__(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: [i for i in sorted(f)]
            [('a', 1), ('c', 3)]

            sage: s = SymmetricFunctions(QQ).schur()                                    # needs sage.combinat
            sage: a = s([2,1]) + s([3])                                                 # needs sage.combinat
            sage: [i for i in sorted(a)]                                                # needs sage.combinat
            [([2, 1], 1), ([3], 1)]
        """
        return iter(self._monomial_coefficients.items())

    def __contains__(self, x):
        """
        Return whether or not a combinatorial object ``x`` indexing a basis
        element is in the support of ``self``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: 'a' in f
            doctest:warning...
            DeprecationWarning: using 'index in vector' is deprecated; use 'index in vector.support()' instead
            See https://github.com/sagemath/sage/issues/34509 for details.
            True
            sage: 'b' in f
            False

            sage: # needs sage.combinat
            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1]) + s([3])
            sage: Partition([2,1]) in a
            True
            sage: Partition([1,1,1]) in a
            False
        """
        deprecation(34509, "using 'index in vector' is deprecated; use 'index in vector.support()' instead")
        return x in self.support()

    def __hash__(self):
        """
        Return the hash value for ``self``.

        The result is cached.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: hash(f) == hash(B['a'] + 3*B['c'])
            True
            sage: hash(f) == hash(B['a'] + 4*B['c'])
            False

            sage: # needs sage.combinat
            sage: F = RootSystem(['A',2]).ambient_space()
            sage: f = F.simple_root(0)
            sage: hash(f) == hash(F.simple_root(0))
            True
            sage: hash(f) == hash(F.simple_root(1))
            False

        This uses the recipe that was proposed for frozendicts in
        :pep:`416` (and adds the hash of the parent). This recipe
        relies on the hash function for frozensets which uses tricks
        to mix the hash values of the items in case they are similar.

        .. TODO::

            It would be desirable to make the hash value depend on the
            hash value of the parent. See :issue:`15959`.
        """
        if not self._hash_set:
            self._hash = hash(frozenset(self._monomial_coefficients.items()))
            self._hash_set = True
        return self._hash

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: loads(dumps(F.an_element())) == F.an_element()
            True
        """
        return (_unpickle_element, (self._parent, self._monomial_coefficients))

    def __setstate__(self, state):
        r"""
        For unpickling old ``CombinatorialFreeModuleElement`` classes.

        See :issue:`22632` and ``register_unpickle_override`` below.

        EXAMPLES::

            sage: loads(b'x\x9c\x95R\xcbn\x131\x14\xd5\x00\r\x89KK\xcb\xa3'
            ....: b'\xbc\xa1\xbc\xd3\xcd,\xe0\x0f\n\xad\xc4\xa2Y\x0c\xb2XZ'
            ....: b'\x8e\xe7N\xe6\x8a\xb1\xa7\xd7\x0f\x91,F\x82E&\xe2\xafq3'
            ....: b'\x13\xa4"X\xb0\xb1}\xae}\xce=\xf7\xc8\xdf\xaf(\'g\x90:o'
            ....: b'\x83\xf2\xc1B\x9a/\x8c\xd4\xa8\x84\xaa\xa4s\xec2\xa2d'
            ....: b'\xcc\xdf\x7f\xa8\xf5\x14\x8d\xf4\xb5EY\x9dZ\x80\xb3:'
            ....: b'\x0f\x15\x88o\xe8K\xa1\xa4\x87Ym\x17)T\xa0\xc1\xf8\x8eH}'
            ....: b'\x17\xd5S\xd3"\xd2\x84^\xf3\xd8?\xf4N:\x01FW\x95\x10\xd3'
            ....: b'\x80\x95G#\x04\x9b\x81\x97\xde[F\xd7:I\x8dN\xad\x17\xa6dU'
            ....: b'\t\r\xbe\xacsF[\xe5\xd6\x9f\x83\x05\x83\x14@X8\xb7\xe0'
            ....: b'\xa2\xb2\xf4X\x1b\x16\x8c\x85<(`4\xe8=v\x13 \xb8\xb43'
            ....: b'\xe8\xd8Y\xbf\xd3\xf5\xee\x89E3s)\x9a\xf8\x10\xac\xb8@'
            ....: b'\xecS\x07\xb2\x8b3\r\x8f2\x1a-\x1bb|\x98\xa3;\x97^\x95'
            ....: b'\xb4\xfd\xd3\xad\xe8FF;|\xbbKJ\xce\xb1\xd6\xb4\xcbG_":'
            ....: b'\x96\x0e\x1d\xdd\\e\xb4W\xee\xf2\xfdS4\xe8\xe1#\xc6\x00'
            ....: b'\\4)+\xda\x8fW\xb7\xf8\xce\xe5To\xb7\x19\xddi\xe9\xeed2'
            ....: b'\xf1\x19\x1d\x1c\xfd\xa0{\xe5\xe0\xff\x93ft\xbf\x1cm\x88'
            ....: b'\x0e\xbcK\x8bu\x7f\x01&h\xb01\x8f\\\xc42\xeb\\\x9d\xfc.~'
            ....: b'\x8e5z\xc0\x939O\x16-=\\6+z\x94\xd1\xe3\xb6\xa1\'c>\xdc'
            ....: b'\xfc\x04zZ\xee\xf1A\xcc\xbc\xc09=\xe3\xc9qX\xd1aF\xcf'
            ....: b'\x1bz\xc1\x0f\xa23S\xeb\xe8F\xa8\x1a\x8a\x02\x15\xc6\xe9'
            ....: b'\x1c\xbdl\xe8\xd58\xaa\xfe%n\xa6\xe5W\x10\x1b@\xafy\xf2n'
            ....: b'\x99\xd1\x9b\xe8\xa2\xec\xcfo\x83k\xa7\xe9/\xc1\xe1\t\x17')
            2*B['x'] + 2*B['y']
        """
        self._set_parent(state[0])
        for k, v in state[1].items():
            setattr(self, k, v)

    def __copy__(self):
        r"""
        Return ``self`` since ``self`` is immutable.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: x = F.an_element()
            sage: copy(x) is x
            True
        """
        return self

    def __deepcopy__(self, memo=None):
        r"""
        Return ``self`` since ``self`` is immutable.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: x = F.an_element()
            sage: deepcopy(x) is x
            True
        """
        return self

    cpdef dict monomial_coefficients(self, bint copy=True):
        """
        Return the internal dictionary which has the combinatorial objects
        indexing the basis as keys and their corresponding coefficients as
        values.

        INPUT:

        - ``copy`` -- boolean (default: ``True``); if ``self`` is internally
          represented by a dictionary ``d``, then make a copy of ``d``.
          If ``False``, then this can cause undesired behavior by
          mutating ``d``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: d = f.monomial_coefficients()
            sage: d['a']
            1
            sage: d['c']
            3

        To run through the monomials of an element, it is better to
        use the idiom::

            sage: for (t,c) in f:
            ....:     print("{} {}".format(t,c))
            a 1
            c 3

        ::

            sage: # needs sage.combinat
            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1])+2*s([3,2])
            sage: d = a.monomial_coefficients()
            sage: type(d)
            <... 'dict'>
            sage: d[ Partition([2,1]) ]
            1
            sage: d[ Partition([3,2]) ]
            2
        """
        if copy:
            return dict(self._monomial_coefficients)
        return self._monomial_coefficients

    def _sorted_items_for_printing(self):
        """
        Return the items (i.e. terms) of ``self``, sorted for printing.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 2*B['c'] + 3 * B['b']
            sage: f._sorted_items_for_printing()
            [('a', 1), ('b', 3), ('c', 2)]
            sage: F.print_options(sorting_reverse=True)
            sage: f._sorted_items_for_printing()
            [('c', 2), ('b', 3), ('a', 1)]
            sage: F.print_options(sorting_reverse=False) #reset to original state

        .. SEEALSO:: :meth:`_repr_`, :meth:`_latex_`, :meth:`print_options`
        """
        print_options = self._parent.print_options()
        v = list(self._monomial_coefficients.items())
        try:
            v.sort(key=lambda monomial_coeff:
                        print_options['sorting_key'](monomial_coeff[0]),
                   reverse=print_options['sorting_reverse'])
        except Exception: # Sorting the output is a plus, but if we can't, no big deal
            pass
        return v

    def _repr_(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='F')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            F['a'] + 2*F['b']
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            ['a'] + 2*['b']
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='', scalar_mult=' ', bracket=False)
            sage: e = F.basis()
            sage: e['a'] + 2*e['b'] # indirect doctest
            'a' + 2 'b'

        Controling the order of terms by providing a comparison
        function on elements of the support::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'],
            ....:                             sorting_reverse=True)
            sage: e = F.basis()
            sage: e['a'] + 3*e['b'] + 2*e['c']
            2*B['c'] + 3*B['b'] + B['a']

            sage: F = CombinatorialFreeModule(QQ, ['ac', 'ba', 'cb'],
            ....:                             sorting_key=lambda x: x[1])
            sage: e = F.basis()
            sage: e['ac'] + 3*e['ba'] + 2*e['cb']
            3*B['ba'] + 2*B['cb'] + B['ac']
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult=self._parent._print_options['scalar_mult'],
                            repr_monomial = self._parent._repr_term,
                            strip_one = True)

    def _ascii_art_(self):
        r"""
        TESTS::

            sage: # needs sage.combinat
            sage: M = QuasiSymmetricFunctions(QQ).M()
            sage: ascii_art(M[1,3]**2)  # indirect doctest
            4*M      + 2*M       + 2*M      + 2*M       + 2*M       + M
                 ***      ******        ***         ***         ***     ******
               ***        *             *        ****         ***      **
               *          *           ***        *           **
               *                      *
            sage: ascii_art(M.zero())
            0
            sage: DA = DescentAlgebra(QQ, 4)
            sage: ascii_art(DA.an_element())
            2*B  + 2*B   + 3*B
               *      **       *
               *      *       **
               *      *       *
               *

        We can get the ascii art when there is no ``one_basis`` method
        (and the basis keys do not compare with ``None``)::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: E2 = L.exterior_power(2)
            sage: ascii_art(E2.an_element())
            2*()/\(5,6,7) + 2*()/\(5,7,6) + 3*()/\(1,2)(3,4)

        We can also get the ascii art when ``one_basis``
        is ``NotImplemented``::

            sage: TL = TemperleyLiebAlgebra(8, -1, QQ)
            sage: C = TL.cellular_basis()
            sage: ascii_art(list(C.basis())[0])
            C
             (     .-. .-. .-. .-.   .-. .-. .-. .-. )
             ( 0,  o o o o o o o o,  o o o o o o o o )
        """
        from sage.misc.repr import coeff_repr
        terms = self._sorted_items_for_printing()
        scalar_mult = self._parent._print_options['scalar_mult']
        repr_monomial = self._parent._ascii_art_term
        strip_one = True

        if repr_monomial is None:
            repr_monomial = str

        chunks = []
        first = True

        if scalar_mult is None:
            scalar_mult = "*"

        one_basis = None
        try:
            if self._parent.one_basis is not NotImplemented:
                one_basis = self._parent.one_basis()
        except (AttributeError, NotImplementedError, ValueError, TypeError):
            pass

        for monomial, c in terms:
            b = repr_monomial(monomial)  # PCR
            if c != 0:
                break_points = []
                coeff = coeff_repr(c, False)
                if coeff != "0":
                    if coeff == "1":
                        coeff = ""
                    elif coeff == "-1":
                        coeff = "-"
                    elif b._l > 0:
                        if one_basis is not None and len(coeff) > 0 and monomial == one_basis and strip_one:
                            b = empty_ascii_art  # ""
                        else:
                            b = AsciiArt([scalar_mult]) + b
                    if not first:
                        if len(coeff) > 0 and coeff[0] == "-":
                            coeff = " - %s" % coeff[1:]
                        else:
                            coeff = " + %s" % coeff
                        break_points = [2]
                    else:
                        coeff = "%s" % coeff
                if coeff:
                    chunks.append(AsciiArt([coeff], break_points))
                if b._l:
                    chunks.append(b)
                first = False
        s = ascii_art(*chunks)
        if first:
            return AsciiArt(["0"])
        if s == empty_ascii_art:
            return AsciiArt(["1"])
        return s

    def _unicode_art_(self):
        r"""
        TESTS::

            sage: M = QuasiSymmetricFunctions(QQ).M()                                   # needs sage.combinat
            sage: unicode_art(M[1,1]**2)  # indirect doctest                            # needs sage.combinat
            6*M   + 2*M    + 2*M    + 2*M    + M
               ┌┐      ┌┬┐       ┌┐       ┌┐     ┌┬┐
               ├┤      ├┼┘      ┌┼┤       ├┤    ┌┼┼┘
               ├┤      ├┤       ├┼┘      ┌┼┤    └┴┘
               ├┤      └┘       └┘       └┴┘
               └┘

        The following test failed before :issue:`26850`::

            sage: unicode_art([M.zero()])  # indirect doctest                           # needs sage.combinat
            [ 0 ]

        We can get the unicode art when there is no ``one_basis`` method
        (and the basis keys do not compare with ``None``)::

            sage: DC3 = groups.permutation.DiCyclic(3)
            sage: L = DC3.regular_representation(QQ, side='left')
            sage: E2 = L.exterior_power(2)
            sage: unicode_art(E2.an_element())
            2*()∧(5,6,7) + 2*()∧(5,7,6) + 3*()∧(1,2)(3,4)

        We can also get the unicode art when ``one_basis``
        is ``NotImplemented``::

            sage: TL = TemperleyLiebAlgebra(8, -1, QQ)
            sage: C = TL.cellular_basis()
            sage: unicode_art(list(C.basis())[0])
            C
             ⎛     ╭─╮ ╭─╮ ╭─╮ ╭─╮   ╭─╮ ╭─╮ ╭─╮ ╭─╮ ⎞
             ⎝ 0,  ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬,  ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⎠
        """
        from sage.misc.repr import coeff_repr
        terms = self._sorted_items_for_printing()
        scalar_mult = self._parent._print_options['scalar_mult']
        repr_monomial = self._parent._unicode_art_term
        strip_one = True

        if repr_monomial is None:
            repr_monomial = str

        chunks = []
        first = True

        if scalar_mult is None:
            scalar_mult = "*"

        one_basis = None
        try:
            if self._parent.one_basis is not NotImplemented:
                one_basis = self._parent.one_basis()
        except (AttributeError, NotImplementedError, ValueError, TypeError):
            pass

        for (monomial, c) in terms:
            b = repr_monomial(monomial)  # PCR
            if c != 0:
                break_points = []
                coeff = coeff_repr(c, False)
                if coeff != "0":
                    if coeff == "1":
                        coeff = ""
                    elif coeff == "-1":
                        coeff = "-"
                    elif b._l > 0:
                        if one_basis is not None and len(coeff) > 0 and monomial == one_basis and strip_one:
                            b = empty_unicode_art  # ""
                        else:
                            b = UnicodeArt([scalar_mult]) + b
                    if not first:
                        if len(coeff) > 0 and coeff[0] == "-":
                            coeff = " - %s" % coeff[1:]
                        else:
                            coeff = " + %s" % coeff
                        break_points = [2]
                    else:
                        coeff = "%s" % coeff
                if coeff:
                    chunks.append(UnicodeArt([coeff], break_points))
                if b._l:
                    chunks.append(b)
                first = False
        s = unicode_art(*chunks)
        if first:
            return UnicodeArt(["0"])
        elif s == empty_unicode_art:
            return UnicodeArt(["1"])
        else:
            return s

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: latex(f)
            B_{a} + 3 B_{c}

        ::

            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)                                    # needs sage.combinat
            sage: a = 2 + QS3([2,1,3])                                                  # needs sage.combinat
            sage: latex(a)  #indirect doctest                                           # needs sage.combinat
            2 [1, 2, 3] + [2, 1, 3]

       ::

            sage: F = CombinatorialFreeModule(QQ, ['a','b'], prefix='beta', latex_prefix='\\beta')
            sage: x = F.an_element()
            sage: x
            2*beta['a'] + 2*beta['b']
            sage: latex(x)
            2 \beta_{a} + 2 \beta_{b}

        Controling the order of terms by providing a comparison
        function on elements of the support::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'],
            ....:                             sorting_reverse=True)
            sage: e = F.basis()
            sage: latex(e['a'] + 3*e['b'] + 2*e['c'])
            2 B_{c} + 3 B_{b} + B_{a}

            sage: F = CombinatorialFreeModule(QQ, ['ac', 'ba', 'cb'],
            ....:                             sorting_key=lambda x: x[1])
            sage: e = F.basis()
            sage: latex(e['ac'] + 3*e['ba'] + 2*e['cb'])
            3 B_{ba} + 2 B_{cb} + B_{ac}
        """
        return repr_lincomb(self._sorted_items_for_printing(),
                            scalar_mult       = self._parent._print_options['scalar_mult'],
                            latex_scalar_mult = self._parent._print_options['latex_scalar_mult'],
                            repr_monomial = self._parent._latex_term,
                            is_latex=True, strip_one=True)

    cpdef _richcmp_(self, other, int op):
        """
        Rich comparison for equal parents.

        EXAMPLES::

            sage: F1 = CombinatorialFreeModule(QQ, [1, 2, 3])
            sage: F2 = CombinatorialFreeModule(QQ, [1, 2, 3], prefix = "g")
            sage: F1.zero() == F1.zero()
            True
            sage: F1.zero() == F1.an_element()
            False
            sage: F1.an_element() == F1.an_element()
            True
            sage: F1.an_element() is None
            False

        ::

            sage: F3 = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: F3.an_element() != F3.an_element()
            False
            sage: F3.an_element() != F3.zero()
            True

        ::

            sage: # needs sage.combinat
            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s([2,1])
            sage: b = s([1,1,1])
            sage: a == b
            False

        .. TODO::

            Currently, if ``self`` and ``other`` do not have the same parent,
            seemingly equal elements do not evaluate equal, since conversions
            between different modules have not been established.

        ::

            sage: F1.zero() == 0
            True
            sage: F1(0)
            0

        ::

            sage: F1.zero() == F2.zero()
            False
            sage: F1(F2.zero())
            Traceback (most recent call last):
            ...
            TypeError: do not know how to make x (= 0) an element of self (=Free module generated by {1, 2, 3} over Rational Field)

            sage: # needs sage.combinat
            sage: F = AlgebrasWithBasis(QQ).example()
            sage: F.one() == 1
            True
            sage: 1 == F.one()
            True
            sage: 2 * F.one() == int(2)
            True
            sage: int(2) == 2 * F.one()
            True

            sage: S = SymmetricFunctions(QQ); s = S.s(); p = S.p()                      # needs sage.combinat
            sage: p[2] == s[2] - s[1, 1]                                                # needs sage.combinat
            True
            sage: p[2] == s[2]                                                          # needs sage.combinat
            False

        This feature is disputable, in particular since it can make
        equality testing costly. It may be removed at some point.

        Equality testing can be a bit tricky when the order of terms
        can vary because their indices are incomparable with
        ``cmp``. The following test did fail before :issue:`12489` ::

            sage: # needs sage.combinat
            sage: F = CombinatorialFreeModule(QQ, Subsets([1,2,3]))
            sage: x = F.an_element()
            sage: (x+F.zero()).terms()  # random
            [2*B[{1}], 3*B[{2}], B[{}]]
            sage: x.terms()             # random
            [2*B[{1}], B[{}], 3*B[{2}]]
            sage: x+F.zero() == x
            True

        TESTS::

            sage: TestSuite(F1).run()
            sage: TestSuite(F).run()                                                    # needs sage.combinat
        """
        cdef IndexedFreeModuleElement elt = <IndexedFreeModuleElement> other

        if self._monomial_coefficients == elt._monomial_coefficients:
            return rich_to_bool(op, 0)

        # Not equal
        if op == Py_EQ:
            return False
        if op == Py_NE:
            return True

        v = sorted(self._monomial_coefficients.items())
        w = sorted(elt._monomial_coefficients.items())
        return richcmp(v, w, op)

    cpdef _add_(self, other):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a'] + 3*B['c']
            B['a'] + 3*B['c']

        ::

            sage: # needs sage.combinat
            sage: s = SymmetricFunctions(QQ).schur()
            sage: s([2,1]) + s([5,4]) # indirect doctest
            s[2, 1] + s[5, 4]
            sage: a = s([2,1]) + 0
            sage: len(a.monomial_coefficients())
            1
        """
        return type(self)(self._parent,
                          add(self._monomial_coefficients,
                              (<IndexedFreeModuleElement>other)._monomial_coefficients))

    cpdef _neg_(self):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] + 3*B['c']
            sage: -f
            -B['a'] - 3*B['c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()                                    # needs sage.combinat
            sage: -s([2,1])  # indirect doctest                                         # needs sage.combinat
            -s[2, 1]
        """
        return type(self)(self._parent, negate(self._monomial_coefficients))

    cpdef _sub_(self, other):
        """
        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a'] - 3*B['c']
            B['a'] - 3*B['c']

        ::

            sage: s = SymmetricFunctions(QQ).schur()                                    # needs sage.combinat
            sage: s([2,1]) - s([5,4])  # indirect doctest                               # needs sage.combinat
            s[2, 1] - s[5, 4]
        """
        return type(self)(self._parent,
                          axpy(-1,
                               (<IndexedFreeModuleElement>other)._monomial_coefficients,
                               self._monomial_coefficients))

    def __getitem__(self, m):
        """
        Return the coefficient of ``m`` in ``self``.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: p = Partition([2,1])
            sage: q = Partition([1,1,1])
            sage: s = SymmetricFunctions(QQ).schur()
            sage: a = s(p)
            sage: a[p]
            1
            sage: a[q]
            0
            sage: a[[2,1]]
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'
        """
        res = self._monomial_coefficients.get(m)
        if res is None:
            return self.base_ring().zero()
        return res

    def _vector_(self, new_base_ring=None, order=None, sparse=False):
        r"""
        Return ``self`` as a vector.

        INPUT:

        - ``new_base_ring`` -- a ring (default: ``None``)
        - ``order`` -- (optional) an ordering of the support of ``self``
        - ``sparse`` -- boolean (default: ``False``); whether to return a sparse
          vector or a dense vector

        OUTPUT: a :func:`FreeModule` vector

        .. WARNING:: This will crash/run forever if ``self`` is infinite dimensional!

        .. SEEALSO::

            - :func:`vector`
            - :meth:`CombinatorialFreeModule.get_order`
            - :meth:`CombinatorialFreeModule.from_vector`
            - :meth:`CombinatorialFreeModule._dense_free_module`

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: f = B['a'] - 3*B['c']
            sage: f._vector_()
            (1, 0, -3)

        One can use equivalently::

            sage: f.to_vector()
            (1, 0, -3)
            sage: vector(f)
            (1, 0, -3)

        More examples::

            sage: # needs sage.combinat
            sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            sage: a = 2*QS3([1,2,3]) + 4*QS3([3,2,1])
            sage: a._vector_()
            (2, 0, 0, 0, 0, 4)
            sage: a.to_vector()
            (2, 0, 0, 0, 0, 4)
            sage: vector(a)
            (2, 0, 0, 0, 0, 4)
            sage: a == QS3.from_vector(a.to_vector())
            True
            sage: a.to_vector(sparse=True)
            (2, 0, 0, 0, 0, 4)

        If ``new_base_ring`` is specified, then a vector over
        ``new_base_ring`` is returned::

            sage: a._vector_(RDF)                                                       # needs sage.combinat
            (2.0, 0.0, 0.0, 0.0, 0.0, 4.0)

        .. NOTE::

            :issue:`13406`: the current implementation has been optimized, at
            the price of breaking the encapsulation for FreeModule
            elements creation, with the following use case as metric,
            on a 2008' Macbook Pro::

                sage: F = CombinatorialFreeModule(QQ, range(10))
                sage: f = F.an_element()
                sage: %timeit f._vector_()   # not tested
                625 loops, best of 3: 17.5 micros per loop

             Other use cases may call for different or further
             optimizations.
        """
        free_module = self._parent._dense_free_module(new_base_ring)
        if sparse:
            free_module = free_module.sparse_module()
        d = self._monomial_coefficients
        zero = free_module.base_ring().zero()
        if sparse:
            if order is None:
                order = {k: i for i,k in enumerate(self._parent.get_order())}
            return free_module.element_class(free_module,
                                             {order[k]: c for k, c in d.items()},
                                             coerce=True, copy=False)
        else:
            if order is None:
                order = self._parent.get_order()
            return free_module.element_class(free_module,
                                             [d.get(m, zero) for m in order],
                                             coerce=True, copy=False)

    to_vector = _vector_

    cpdef _acted_upon_(self, scalar, bint self_on_left):
        """
        Return the action of ``scalar`` (an element of the base ring) on
        ``self``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a','b','c'])
            sage: B = F.basis()
            sage: B['a']*(1/2)  # indirect doctest
            1/2*B['a']
            sage: B['a']/2
            1/2*B['a']
            sage: B['a']*2      # indirect doctest
            2*B['a']
            sage: B['a']*int(2) # indirect doctest
            2*B['a']

            sage: 1/2*B['a']
            1/2*B['a']
            sage: 2*B['a']      # indirect doctest
            2*B['a']
            sage: int(2)*B['a'] # indirect doctest
            2*B['a']

        TESTS::

            sage: coercion_model.get_action(F, QQ, operator.mul)
            Right scalar multiplication by Rational Field on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: coercion_model.get_action(QQ, F, operator.mul)
            Left scalar multiplication by Rational Field on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: coercion_model.get_action(F, ZZ, operator.mul)
            Right scalar multiplication by Integer Ring on Free module generated by {'a', 'b', 'c'} over Rational Field
            sage: print(coercion_model.get_action(F, F, operator.mul))
            None

        This also works when a coercion of the coefficient is needed, for
        example with polynomials or fraction fields (:issue:`8832`)::

            sage: P.<q> = QQ['q']
            sage: V = CombinatorialFreeModule(P, Permutations())
            sage: el = V(Permutation([3,1,2]))
            sage: (3/2)*el
            3/2*B[[3, 1, 2]]

            sage: # needs sage.combinat
            sage: P.<q> = QQ['q']
            sage: F = FractionField(P)
            sage: V = CombinatorialFreeModule(F, Words())
            sage: w = Words()('abc')
            sage: (1+q)*V(w)
            (q+1)*B[word: abc]
            sage: ((1+q)/q)*V(w)
            ((q+1)/q)*B[word: abc]

        .. TODO::

            Add non commutative tests.
        """
        # With the current design, the coercion model does not have
        # enough information to detect a priori that this method only
        # accepts scalars; so it tries on some elements(), and we need
        # to make sure to report an error.
        if isinstance(scalar, Element) and parent(scalar) is not self.base_ring():
            # Temporary needed by coercion (see Polynomial/FractionField tests).
            if self.base_ring().has_coerce_map_from(parent(scalar)):
                scalar = self.base_ring()( scalar )
            else:
                return None

        return type(self)(self._parent,
                          scal(scalar, self._monomial_coefficients,
                               factor_on_left=not self_on_left))

    cpdef _lmul_(self, Element right):
        """
        For backward compatibility.

        EXAMPLES::

            sage: C = CombinatorialFreeModule(QQ, [1,2,3])
            sage: C.an_element()._lmul_(2)
            4*B[1] + 4*B[2] + 6*B[3]
        """
        return self._acted_upon_(right, True)

    cpdef _rmul_(self, Element left):
        """
        For backward compatibility.

        EXAMPLES::

            sage: C = CombinatorialFreeModule(QQ, [1,2,3])
            sage: C.an_element()._rmul_(2)
            4*B[1] + 4*B[2] + 6*B[3]
        """
        return self._acted_upon_(left, False)

    def __truediv__(left, x):
        """
        Division by coefficients.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, [1,2,3])
            sage: x = F._from_dict({1:2, 2:3})
            sage: from operator import truediv
            sage: truediv(x, 2)
            B[1] + 3/2*B[2]

        ::

            sage: F = CombinatorialFreeModule(QQ, [1,2,3])
            sage: B = F.basis()
            sage: f = 2*B[2] + 4*B[3]
            sage: truediv(f, 2)
            B[2] + 2*B[3]

        TESTS::

            sage: truediv(x, x)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 2*B[1] + 3*B[2] to a rational
            sage: truediv("hello", x)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for /: 'str' and 'CombinatorialFreeModule_with_category.element_class'

            sage: L = LazyPowerSeriesRing(QQ, 't')
            sage: t = L.gen()
            sage: F = algebras.Free(L, ['A', 'B'])
            sage: A, B = F.gens()
            sage: f = t*A + t**2*B/2
        """
        if not isinstance(left, IndexedFreeModuleElement):
            return NotImplemented

        cdef IndexedFreeModuleElement self = <IndexedFreeModuleElement>left
        F = self._parent
        B = self.base_ring()
        D = self._monomial_coefficients
        if B not in _Fields:
            return type(self)(F, {k: c._divide_if_possible(x)
                                  for k, c in D.items()})

        x_inv = B(x) ** -1
        return type(self)(F, scal(x_inv, D))

    def _magma_init_(self, magma):
        r"""
        Convert ``self`` to Magma.
        """
        # Get a reference to Magma version of parent.
        R = magma(self.parent()).name()
        # use dict {key: coefficient}.
        return '+'.join(f"({c._magma_init_(magma)})*{R}.{m._magma_init_(magma)}"
                        for m, c in self.monomial_coefficients().items())


def _unpickle_element(C, d):
    """
    Unpickle an element in ``C`` given by ``d``.

    EXAMPLES::

        sage: from sage.modules.with_basis.indexed_element import _unpickle_element
        sage: C = CombinatorialFreeModule(QQ, [1,2,3])
        sage: _unpickle_element(C, {1: -2, 3: -12})
        -2*B[1] - 12*B[3]
    """
    return C._from_dict(d, coerce=False, remove_zeros=False)


# Handle old CombinatorialFreeModuleElement pickles, see Issue #22632
from sage.misc.persist import register_unpickle_override
register_unpickle_override("sage.combinat.free_module",
                           "CombinatorialFreeModuleElement",
                           IndexedFreeModuleElement)
