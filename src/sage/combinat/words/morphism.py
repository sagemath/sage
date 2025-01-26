r"""
Word morphisms/substitutions

This module implements morphisms over finite and infinite words.

AUTHORS:

- Sébastien Labbé (2007-06-01): initial version
- Sébastien Labbé (2008-07-01): merged into sage-words
- Sébastien Labbé (2008-12-17): merged into sage
- Sébastien Labbé (2009-02-03): words next generation
- Sébastien Labbé (2009-11-20): allowing the choice of the
  datatype of the image. Doc improvements.
- Stepan Starosta (2012-11-09): growing letters

EXAMPLES:

Creation of a morphism from a dictionary or a string::

    sage: n = WordMorphism({0:[0,2,2,1],1:[0,2],2:[2,2,1]})

::

    sage: m = WordMorphism('x->xyxsxss,s->xyss,y->ys')

::

    sage: n
    WordMorphism: 0->0221, 1->02, 2->221
    sage: m
    WordMorphism: s->xyss, x->xyxsxss, y->ys

The codomain may be specified::

    sage: WordMorphism({0:[0,2,2,1],1:[0,2],2:[2,2,1]}, codomain=Words([0,1,2,3,4]))
    WordMorphism: 0->0221, 1->02, 2->221

Power of a morphism::

    sage: n^2
    WordMorphism: 0->022122122102, 1->0221221, 2->22122102

Image under a morphism::

    sage: m('y')
    word: ys
    sage: m('xxxsy')
    word: xyxsxssxyxsxssxyxsxssxyssys

Iterated image under a morphism::

    sage: m('y', 3)
    word: ysxyssxyxsxssysxyssxyss

See more examples in the documentation of the call method
(``m.__call__?``).

Infinite fixed point of morphism::

    sage: fix = m.fixed_point('x')
    sage: fix
    word: xyxsxssysxyxsxssxyssxyxsxssxyssxyssysxys...
    sage: fix.length()
    +Infinity

Incidence matrix::

    sage: matrix(m)                                                                     # needs sage.modules
    [2 3 1]
    [1 3 0]
    [1 1 1]

Many other functionalities...::

    sage: m.is_identity()
    False
    sage: m.is_endomorphism()
    True
"""
# ****************************************************************************
#       Copyright (C) 2008 Sebastien Labbe <slabqc@gmail.com>
#                     2018 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Iterable
from sage.misc.callable_dict import CallableDict
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.lazy_list import lazy_list
from sage.sets.set import Set
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer import Integer
from sage.combinat.words.word import FiniteWord_class
from sage.combinat.words.words import FiniteWords, FiniteOrInfiniteWords

lazy_import('sage.modules.free_module_element', 'vector')
lazy_import('sage.matrix.constructor', 'Matrix')


def get_cycles(f, domain):
    r"""
    Return the list of cycles of the function ``f`` contained in ``domain``.

    INPUT:

    - ``f`` -- function

    - ``domain`` -- iterable, a subdomain of the domain of definition of ``f``

    EXAMPLES::

        sage: from sage.combinat.words.morphism import get_cycles
        sage: get_cycles(lambda i: (i+1)%3, [0,1,2])
        [(0, 1, 2)]
        sage: get_cycles(lambda i: [0,0,0][i], [0,1,2])
        [(0,)]
        sage: get_cycles(lambda i: [1,1,1][i], [0,1,2])
        [(1,)]
        sage: get_cycles(lambda i: [2,3,0][i], [0,1,2])
        [(0, 2)]
        sage: d = {'a': 'a', 'b': 'b'}
        sage: get_cycles(d.__getitem__, 'ba')
        [('b',), ('a',)]
    """
    cycles = []
    not_seen = set(domain)
    for a in domain:
        if a not in not_seen:
            continue
        cycle = [a]
        b = f(a)
        not_seen.remove(a)
        while b in not_seen:
            not_seen.remove(b)
            cycle.append(b)
            b = f(b)
        if b in cycle:
            cycles.append(tuple(cycle[cycle.index(b):]))

    return cycles


class PeriodicPointIterator:
    r"""
    (Lazy) constructor of the periodic points of a word morphism.

    This class is mainly used in :class:`WordMorphism.periodic_point` and
    :class:`WordMorphism.periodic_points`.

    EXAMPLES::

        sage: from sage.combinat.words.morphism import PeriodicPointIterator
        sage: s = WordMorphism('a->bacca,b->cba,c->aab')
        sage: p = PeriodicPointIterator(s, ['a','b','c'])
        sage: p._cache[0]
        lazy list ['a', 'a', 'b', ...]
        sage: p._cache[1]
        lazy list ['b', 'a', 'c', ...]
        sage: p._cache[2]
        lazy list ['c', 'b', 'a', ...]
    """
    def __init__(self, m, cycle):
        r"""
        INPUT:

        - ``m`` -- a word morphism

        - ``cycle`` -- a cycle of letters under the morphism

        TESTS::

            sage: from sage.combinat.words.morphism import PeriodicPointIterator
            sage: s = WordMorphism('a->bacca,b->cba,c->aab')
            sage: p = PeriodicPointIterator(s, ['a','b','c'])
            sage: pp = loads(dumps(p))
            sage: pp._cache[0]
            lazy list ['a', 'a', 'b', ...]
        """
        self._m = m            # for pickling only
        self._image = m.image
        self._cycle = tuple(cycle)
        self._cache = [lazy_list(self.get_iterator(i)) for i in range(len(cycle))]

    def __reduce__(self):
        r"""
        TESTS::

            sage: from sage.combinat.words.morphism import PeriodicPointIterator
            sage: s = WordMorphism('a->bacca,b->cba,c->aab')
            sage: p = PeriodicPointIterator(s, ['a','b','c'])
            sage: p.__reduce__()
            (<class 'sage.combinat.words.morphism.PeriodicPointIterator'>,
             (WordMorphism: a->bacca, b->cba, c->aab, ('a', 'b', 'c')))
        """
        return PeriodicPointIterator, (self._m, self._cycle)

    @cached_method
    def get_iterator(self, i):
        r"""
        Internal method.

        EXAMPLES::

            sage: from sage.combinat.words.morphism import PeriodicPointIterator
            sage: s = WordMorphism('a->bacca,b->cba,c->aab')
            sage: p = PeriodicPointIterator(s, ['a','b','c'])
            sage: p.get_iterator(0)
            <generator object ...get_iterator at ...>
        """
        j = (i - 1) % len(self._cycle)
        for a in self._image(self._cycle[j]):
            yield a
        u = iter(self._cache[j])
        next(u)
        while True:
            for a in self._image(next(u)):
                yield a


class WordMorphism(SageObject):
    r"""
    WordMorphism class.

    INPUT:

    - ``data`` -- dictionary or string or an instance of WordMorphism, the map
      giving the image of letters
    - ``domain`` -- (optional:``None``) set of words over a given
      alphabet. If ``None``, the domain alphabet is computed from ``data``
      and is *sorted*.
    - ``codomain`` -- (optional:``None``) set of words over a given
      alphabet. If ``None``, the codomain alphabet is computed from
      ``data`` and is *sorted*.

    .. NOTE::

        When the domain or the codomain are not explicitly given, it is
        expected that the letters are comparable because the alphabets of
        the domain and of the codomain are sorted.

    EXAMPLES:

    From a dictionary::

        sage: n = WordMorphism({0:[0,2,2,1],1:[0,2],2:[2,2,1]})
        sage: n
        WordMorphism: 0->0221, 1->02, 2->221

    From a string with ``'->'`` as separation::

        sage: m = WordMorphism('x->xyxsxss,s->xyss,y->ys')
        sage: m
        WordMorphism: s->xyss, x->xyxsxss, y->ys
        sage: m.domain()
        Finite words over {'s', 'x', 'y'}
        sage: m.codomain()
        Finite words over {'s', 'x', 'y'}

    Specifying the domain and codomain::

        sage: W = FiniteWords([0,1,2])
        sage: d = {0:[0,1], 1:[0,1,0], 2:[0]}
        sage: m = WordMorphism(d, domain=W, codomain=W)
        sage: m([0]).parent()
        Finite words over {0, 1, 2}

    When the alphabet is non-sortable, the domain and/or codomain must be
    explicitly given::

        sage: W = FiniteWords(['a',6])
        sage: d = {'a':['a',6,'a'],6:[6,6,6,'a']}
        sage: WordMorphism(d, domain=W, codomain=W)
        WordMorphism: 6->666a, a->a6a

    TESTS::

        sage: wm = WordMorphism('a->ab,b->ba')
        sage: wm == loads(dumps(wm))
        True
    """
    def __init__(self, data, domain=None, codomain=None):
        r"""
        Construction of the morphism.

        EXAMPLES:

        1. If data is a str::

            sage: WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: WordMorphism('a->abc,b->bca,c->cab')
            WordMorphism: a->abc, b->bca, c->cab
            sage: WordMorphism('a->abdsf,b->hahdad,c->asdhasd')
            WordMorphism: a->abdsf, b->hahdad, c->asdhasd
            sage: WordMorphism('(->(),)->)(')
            WordMorphism: (->(), )->)(
            sage: WordMorphism('a->53k,b->y5?,$->49i')
            WordMorphism: $->49i, a->53k, b->y5?

        An erasing morphism::

            sage: WordMorphism('a->ab,b->')
            WordMorphism: a->ab, b->

        Use the arrows ('->') correctly::

            sage: WordMorphism('a->ab,b-')
            Traceback (most recent call last):
            ...
            ValueError: the second and third characters must be '->' (not '-')
            sage: WordMorphism('a->ab,b')
            Traceback (most recent call last):
            ...
            ValueError: the second and third characters must be '->' (not '')
            sage: WordMorphism('a->ab,a-]asdfa')
            Traceback (most recent call last):
            ...
            ValueError: the second and third characters must be '->' (not '-]')

        Each letter must be defined only once::

            sage: WordMorphism('a->ab,a->ba')
            Traceback (most recent call last):
            ...
            ValueError: the image of 'a' is defined twice

        2. From a dictionary::

            sage: WordMorphism({"a":"ab","b":"ba"})
            WordMorphism: a->ab, b->ba
            sage: WordMorphism({2:[4,5,6],3:[1,2,3]})
            WordMorphism: 2->456, 3->123

        The image of a letter can be a set, but the order is not
        preserved::

            sage: WordMorphism({2:[4,5,6],3:set([4,1,8])})  # random results
            WordMorphism: 2->456, 3->814

        If the image of a letter is not iterable, it is considered as a
        letter::

            sage: WordMorphism({0:1, 1:0})
            WordMorphism: 0->1, 1->0
            sage: WordMorphism({0:123, 1:789})
            WordMorphism: 0->123, 1->789
            sage: WordMorphism({2:[4,5,6], 3:123})
            WordMorphism: 2->456, 3->123

        3. From a WordMorphism::

            sage: WordMorphism(WordMorphism('a->ab,b->ba'))
            WordMorphism: a->ab, b->ba

        TESTS::

            sage: WordMorphism(',,,a->ab,,,b->ba,,')
            WordMorphism: a->ab, b->ba

            sage: WordMorphism({(1,2):'ab', 'a': ['c', (1,2), 'a']})
            WordMorphism: (1, 2)->ab, a->c,(1, 2),a

            sage: WordMorphism({'a':'a'}, domain=FiniteWords('ab'))
            Traceback (most recent call last):
            ...
            ValueError: invalid input; the keys of the dictionary must coincide with the domain alphabet
            sage: WordMorphism({'a':'a', 'b':'b'}, domain=FiniteWords('a'))
            Traceback (most recent call last):
            ...
            ValueError: invalid input; the keys of the dictionary must coincide with the domain alphabet
        """
        if isinstance(data, WordMorphism):
            self._domain = data._domain
            self._codomain = data._codomain
            self._morph = data._morph
        else:
            if isinstance(data, str):
                data = self._build_dict(data)
            elif not isinstance(data, dict):
                raise NotImplementedError

            if codomain is None:
                codomain = self._build_codomain(data)

            if isinstance(codomain, FiniteOrInfiniteWords):
                codomain = codomain.finite_words()
            elif not isinstance(codomain, FiniteWords):
                raise TypeError("the codomain must be a set of finite words")
            self._codomain = codomain

            self._morph = {}

            dom_alph = []
            for key, val in data.items():
                dom_alph.append(key)
                if val in codomain.alphabet():
                    self._morph[key] = codomain([val])
                else:
                    self._morph[key] = codomain(val)

            if domain is not None:
                if isinstance(domain, FiniteOrInfiniteWords):
                    domain = domain.finite_words()
                elif not isinstance(domain, FiniteWords):
                    raise TypeError("the domain must be a set of finite words")
                A = domain.alphabet()
                if len(self._morph) != A.cardinality() or not all(a in A for a in self._morph):
                    raise ValueError('invalid input; the keys of the dictionary must coincide with the domain alphabet')
            else:
                try:
                    dom_alph.sort()
                except TypeError:
                    dom_alph.sort(key=str)
                domain = FiniteWords(dom_alph)
            self._domain = domain

    def _build_dict(self, s):
        r"""
        Parse the string input to WordMorphism and build the dictionary
        it represents.

        TESTS::

            sage: wm = WordMorphism('a->ab,b->ba')
            sage: wm._build_dict('a->ab,b->ba') == {'a': 'ab', 'b': 'ba'}
            True
            sage: wm._build_dict('a->ab,a->ba')
            Traceback (most recent call last):
            ...
            ValueError: the image of 'a' is defined twice
            sage: wm._build_dict('a->ab,b>ba')
            Traceback (most recent call last):
            ...
            ValueError: the second and third characters must be '->' (not '>b')
        """
        tmp_dict = {}
        for fleche in s.split(','):
            if len(fleche) == 0:
                continue

            if len(fleche) < 3 or fleche[1:3] != '->':
                raise ValueError("the second and third characters must be '->' (not '%s')" % fleche[1:3])

            lettre = fleche[0]
            image = fleche[3:]

            if lettre in tmp_dict:
                raise ValueError("the image of %r is defined twice" % lettre)

            tmp_dict[lettre] = image
        return tmp_dict

    def _build_codomain(self, data):
        r"""
        Return a Words domain containing all the letters in the values of
        data (which must be a dictionary).

        TESTS:

        If the image of all the letters are iterable::

            sage: wm = WordMorphism('a->ab,b->ba')
            sage: wm._build_codomain({'a': 'ab', 'b': 'ba'})
            Finite words over {'a', 'b'}
            sage: wm._build_codomain({'a': 'dcb', 'b': 'a'})
            Finite words over {'a', 'b', 'c', 'd'}
            sage: wm._build_codomain({2:[4,5,6],3:[1,2,3]})
            Finite words over {1, 2, 3, 4, 5, 6}
            sage: wm._build_codomain({2:[4,5,6],3:set([4,1,8])})
            Finite words over {1, 4, 5, 6, 8}

        If the image of a letter is not iterable, it is considered as
        a letter::

            sage: wm._build_codomain({2:[4,5,6],3:123})
            Finite words over {4, 5, 6, 123}
            sage: wm._build_codomain({0:1, 1:0, 2:2})
            Finite words over {0, 1, 2}
        """
        codom_alphabet = set()
        for val in data.values():
            try:
                it = iter(val)
            except TypeError:
                it = [val]
            codom_alphabet.update(it)
        try:
            codom_alphabet = sorted(codom_alphabet)
        except TypeError:
            codom_alphabet = sorted(codom_alphabet, key=str)
        return FiniteWords(codom_alphabet)

    @cached_method
    def __hash__(self):
        r"""
        TESTS::

            sage: hash(WordMorphism('a->ab,b->ba')) # random
            7211091143079804375
        """
        return hash(tuple((k, v) for k, v in self._morph.items())) ^ hash(self._codomain)

    def __eq__(self, other):
        r"""
        Return ``True`` if ``self`` is equal to ``other``.

        EXAMPLES::

            sage: n = WordMorphism('a->a,b->aa,c->aaa')
            sage: n**3 == n**1
            True
            sage: WordMorphism('b->ba,a->ab') == WordMorphism('a->ab,b->ba')
            True
            sage: WordMorphism('b->ba,a->ab') == WordMorphism({"a":"ab","b":"ba"})
            True
            sage: m = WordMorphism({0:[1,2,3],1:[4,5,6]}); m
            WordMorphism: 0->123, 1->456
            sage: o = WordMorphism('0->123,1->456'); o
            WordMorphism: 0->123, 1->456
            sage: m == o
            False

        TESTS:

        Check that equality depends on the codomain::

            sage: m = WordMorphism('a->a,b->aa,c->aaa')
            sage: n = WordMorphism('a->a,b->aa,c->aaa', codomain=Words('abc'))
            sage: m == n
            False
        """
        if not isinstance(other, WordMorphism):
            return False
        return self._morph == other._morph and self._codomain == other._codomain

    def __ne__(self, other):
        r"""
        Return whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->baba')
            sage: n = WordMorphism('a->ab,b->baba')
            sage: o = WordMorphism('a->ab,b->bab')
            sage: m != n
            False
            sage: n != o
            True

        This solves :issue:`12475`::

            sage: s = WordMorphism('1->121,2->131,3->4,4->1')
            sage: s == s.reversal()
            True
            sage: s != s.reversal()
            False
        """
        return not self == other

    def __repr__(self) -> str:
        r"""
        Return the string representation of the morphism.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->ba')
            WordMorphism: a->ab, b->ba
            sage: WordMorphism({0:[0,1],1:[1,0]})
            WordMorphism: 0->01, 1->10

        TESTS::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: repr(s)
            'WordMorphism: a->ab, b->ba'
        """
        return "WordMorphism: %s" % str(self)

    def __str__(self) -> str:
        r"""
        Return the morphism in str.

        EXAMPLES::

            sage: print(WordMorphism('a->ab,b->ba'))
            a->ab, b->ba
            sage: print(WordMorphism({0:[0,1],1:[1,0]}))
            0->01, 1->10

        The output is sorted to make it unique::

            sage: print(WordMorphism('b->ba,a->ab'))
            a->ab, b->ba

        The str method is used for string formatting::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: "Here is a map : %s" % s
            'Here is a map : a->ab, b->ba'

        ::

            sage: s = WordMorphism({1:[1,2],2:[1]})
            sage: s.dual_map()                                                          # needs sage.modules
            E_1^*(1->12, 2->1)

        TESTS::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: str(s)
            'a->ab, b->ba'
        """
        L = [str(lettre) + '->' + image.string_rep()
             for lettre, image in self._morph.items()]
        return ', '.join(sorted(L))

    def __call__(self, w, order=1):
        r"""
        Return the image of ``w`` under ``self`` to the given order.

        INPUT:

        - ``w`` -- word or sequence in the domain of self

        - ``order`` -- integer or plus ``Infinity`` (default: 1)

        OUTPUT: ``word`` -- ``order``-th iterated image under ``self`` of ``w``

        EXAMPLES:

        The image of a word under a morphism:

        1. The image of a finite word under a morphism::

            sage: tm = WordMorphism ('a->ab,b->ba')
            sage: tm('a')
            word: ab
            sage: tm('aabababb')
            word: ababbaabbaabbaba

        2. The iterated image of a word::

            sage: tm('a', 2)
            word: abba
            sage: tm('aba', 3)
            word: abbabaabbaababbaabbabaab

        3. The infinitely iterated image of a letter::

            sage: tm('a', oo)
            word: abbabaabbaababbabaababbaabbabaabbaababba...

        4. The image of an infinite word::

            sage: t = words.ThueMorseWord()
            sage: n = WordMorphism({0:[0, 1], 1:[1, 0]})
            sage: n(t)
            word: 0110100110010110100101100110100110010110...
            sage: n(t, 3)
            word: 0110100110010110100101100110100110010110...
            sage: n(t)[:1000] == t[:1000]
            True

        The Fibonacci word::

            sage: w = words.FibonacciWord()
            sage: m = WordMorphism({0:'a', 1:'b'})
            sage: m(w)
            word: abaababaabaababaababaabaababaabaababaaba...
            sage: f = words.FibonacciWord('ab')
            sage: f[:1000] == m(w)[:1000]
            True

        ::

            sage: w = words.FibonacciWord("ab")
            sage: m = WordMorphism('a->01,b->101')
            sage: m(w)
            word: 0110101011010110101011010101101011010101...

        The word must be in the domain of self::

            sage: tm('0021')
            Traceback (most recent call last):
            ...
            ValueError: 0 not in alphabet

        The order must be a nonnegative integer or plus Infinity::

            sage: tm('a', -1)
            Traceback (most recent call last):
            ...
            TypeError: order (-1) must be a nonnegative integer or plus Infinity
            sage: tm('a', 6.7)
            Traceback (most recent call last):
            ...
            TypeError: order (6.7...) must be a nonnegative integer or plus Infinity

        Only the first letter is considered for infinitely iterated image of
        a word under a morphism::

            sage: tm('aba',oo)
            word: abbabaabbaababbabaababbaabbabaabbaababba...

        The morphism ``self`` must be prolongable on the given letter for infinitely
        iterated image::

            sage: m = WordMorphism('a->ba,b->ab')
            sage: m('a', oo)
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on a

        The empty word is fixed by any morphism for all natural
        powers::

            sage: phi = WordMorphism('a->ab,b->a')
            sage: phi(Word())
            word:
            sage: phi(Word(), oo)
            word:
            sage: it = iter([])
            sage: phi(it, oo)
            word:

        TESTS::

            sage: for i in range(6):
            ....:   tm('a', i)
            word: a
            word: ab
            word: abba
            word: abbabaab
            word: abbabaabbaababba
            word: abbabaabbaababbabaababbaabbabaab
            sage: m = WordMorphism('a->,b->')
            sage: m('')
            word:

        When the input is a finite word, the output is another
        finite word::

            sage: w = m('aabb')
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_char'>

            sage: w == loads(dumps(w))
            True
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile(suffix='.sobj') as f:
            ....:     save(w, filename=f.name)
        """
        if order == 1:
            D = self.domain()
            C = self.codomain()
            if isinstance(w, (tuple, str, list)):
                w = D(w)

            if isinstance(w, FiniteWord_class):
                im = C()
                for a in w:
                    im += self._morph[a]
                return im

            if isinstance(w, Iterable):
                pass
            elif w in self._domain.alphabet():
                return self._morph[w]
            else:
                raise TypeError("do not know how to handle an input (=%s) that is not iterable or not in the domain alphabet" % w)

            # here we assume (maybe wrongly) that the length is infinite
            parent = self.codomain().shift()
            iterator = (x for y in w for x in self._morph[y])
            parent = parent.shift()
            return parent(iterator)

        elif order is Infinity:
            if isinstance(w, (tuple, str, list, FiniteWord_class)):
                if len(w) == 0:
                    return self.codomain()()
                else:
                    letter = w[0]
            elif isinstance(w, Iterable):
                try:
                    letter = next(w)
                except StopIteration:
                    return self.codomain()()
            elif w in self._domain.alphabet():
                letter = w
            else:
                raise TypeError("do not know how to handle an input (=%s) that is not iterable or not in the domain alphabet" % w)
            return self.fixed_point(letter=letter)

        elif isinstance(order, (int, Integer)) and order > 1:
            return self(self(w, order - 1))

        elif order == 0:
            return self._domain(w)

        else:
            raise TypeError("order (%s) must be a nonnegative integer or plus Infinity" % order)

    def latex_layout(self, layout=None):
        r"""
        Get or set the actual latex layout (oneliner vs array).

        INPUT:

        - ``layout`` -- string (default: ``None``); can take one of the
          following values:

          - ``None`` -- returns the actual latex layout; by default, the
            layout is ``'array'``
          - ``'oneliner'`` -- set the layout to ``'oneliner'``
          - ``'array'`` -- set the layout to ``'array'``

        EXAMPLES::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: s.latex_layout()
            'array'
            sage: s.latex_layout('oneliner')
            sage: s.latex_layout()
            'oneliner'
        """
        if layout is None:
            # return the layout
            if not hasattr(self, '_latex_layout'):
                self._latex_layout = 'array'
            return self._latex_layout
        else:
            # change the layout
            self._latex_layout = layout

    def _latex_(self):
        r"""
        Return the latex representation of the morphism.

        Use :meth:`latex_layout` to change latex layout (oneliner vs
        array). The default is a latex array.

        EXAMPLES::

            sage: s = WordMorphism('a->ab,b->ba')
            sage: s._latex_()
            \begin{array}{l}
            a \mapsto ab\\
            b \mapsto ba
            \end{array}

        Change the latex layout to a one liner::

            sage: s.latex_layout('oneliner')
            sage: s._latex_()
            a \mapsto ab,b \mapsto ba

        TESTS:

        Unknown latex style::

            sage: s.latex_layout('tabular')
            sage: s._latex_()
            Traceback (most recent call last):
            ...
            ValueError: unknown latex_layout(=tabular)
        """
        from sage.misc.latex import LatexExpr
        A = self.domain().alphabet()
        latex_layout = self.latex_layout()
        if latex_layout == 'oneliner':
            lines = (fr"{a} \mapsto {self.image(a)}" for a in A)
            return LatexExpr(r','.join(lines))
        if latex_layout == 'array':
            s = r"\begin{array}{l}" + '\n'
            lines = (fr"{a} \mapsto {self.image(a)}" for a in A)
            s += '\\\\\n'.join(lines)
            s += '\n' + r"\end{array}"
            return LatexExpr(s)
        raise ValueError('unknown latex_layout(=%s)' % latex_layout)

    def __mul__(self, other):
        r"""
        Return the morphism ``self``\*``other``.

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo*m
            WordMorphism: a->aba, b->aab
            sage: fibo*fibo
            WordMorphism: a->aba, b->ab
            sage: m*fibo
            WordMorphism: a->abba, b->ab

        ::

            sage: n = WordMorphism('a->a,b->aa,c->aaa')
            sage: p1 = n*m
            sage: p1
            WordMorphism: a->aaa, b->aaa
            sage: p1.domain()
            Finite words over {'a', 'b'}
            sage: p1.codomain()
            Finite words over {'a'}

        ::

            sage: p2 = m*n
            sage: p2
            WordMorphism: a->ab, b->abab, c->ababab
            sage: p2.domain()
            Finite words over {'a', 'b', 'c'}
            sage: p2.codomain()
            Finite words over {'a', 'b'}

        ::

            sage: m = WordMorphism('0->a,1->b')
            sage: n = WordMorphism('a->c,b->e',codomain=Words('abcde'))
            sage: p = n * m
            sage: p.codomain()
            Finite words over {'a', 'b', 'c', 'd', 'e'}

        TESTS::

            sage: m = WordMorphism('a->b,b->c,c->a')
            sage: WordMorphism('')*m
            Traceback (most recent call last):
            ...
            KeyError: 'b'
            sage: m * WordMorphism('')
            WordMorphism:
        """
        return WordMorphism({key: self(w) for key, w in other._morph.items()},
                            codomain=self.codomain())

    def __pow__(self, exp):
        r"""
        Return the power of ``self`` with exponent = ``exp``.

        INPUT:

        - ``exp`` -- positive integer

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: m^1
            WordMorphism: a->ab, b->ba
            sage: m^2
            WordMorphism: a->abba, b->baab
            sage: m^3
            WordMorphism: a->abbabaab, b->baababba

        The exponent must be a positive integer::

            sage: m^1.5
            Traceback (most recent call last):
            ...
            ValueError: exponent (1.5...) must be an integer
            sage: m^-2
            Traceback (most recent call last):
            ...
            ValueError: exponent (-2) must be strictly positive

        When ``self`` is not an endomorphism::

            sage: n = WordMorphism('a->ba,b->abc')
            sage: n^2
            Traceback (most recent call last):
            ...
            KeyError: 'c'
        """
        # If exp is not an integer
        if not isinstance(exp, (int, Integer)):
            raise ValueError("exponent (%s) must be an integer" % exp)

        # If exp is negative
        elif exp <= 0:
            raise ValueError("exponent (%s) must be strictly positive" % exp)

        # Base of induction
        elif exp == 1:
            return self

        else:
            nexp = int(exp // 2)
            over = exp % 2
            res = (self * self)**nexp
            if over == 1:
                res *= self
            return res

    def extend_by(self, other):
        r"""
        Return ``self`` extended by ``other``.

        Let `\varphi_1:A^*\rightarrow B^*` and `\varphi_2:C^*\rightarrow D^*`
        be two morphisms. A morphism `\mu:(A\cup C)^*\rightarrow (B\cup D)^*`
        corresponds to `\varphi_1` *extended by* `\varphi_2` if
        `\mu(a)=\varphi_1(a)` if `a\in A` and `\mu(a)=\varphi_2(a)` otherwise.

        INPUT:

        - ``other`` -- a WordMorphism

        OUTPUT: WordMorphism

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: n = WordMorphism({'0':'1','1':'0','a':'5'})
            sage: m.extend_by(n)
            WordMorphism: 0->1, 1->0, a->ab, b->ba
            sage: n.extend_by(m)
            WordMorphism: 0->1, 1->0, a->5, b->ba
            sage: m.extend_by(m)
            WordMorphism: a->ab, b->ba

        TESTS::

            sage: m.extend_by(WordMorphism({})) == m
            True
            sage: m.extend_by(WordMorphism('')) == m
            True

        ::

            sage: m.extend_by(4)
            Traceback (most recent call last):
            ...
            TypeError: other (=4) is not a WordMorphism
        """
        if not isinstance(other, WordMorphism):
            raise TypeError("other (=%s) is not a WordMorphism" % other)

        nv = dict(other._morph)
        nv.update(self._morph)
        return WordMorphism(nv)

    def restrict_domain(self, alphabet):
        r"""
        Return a restriction of ``self`` to the given alphabet.

        INPUT:

        - ``alphabet`` -- an iterable

        OUTPUT: WordMorphism

        EXAMPLES::

            sage: m = WordMorphism('a->b,b->a')
            sage: m.restrict_domain('a')
            WordMorphism: a->b
            sage: m.restrict_domain('')
            WordMorphism:
            sage: m.restrict_domain('A')
            WordMorphism:
            sage: m.restrict_domain('Aa')
            WordMorphism: a->b

        The input alphabet must be iterable::

            sage: m.restrict_domain(66)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
        """
        return WordMorphism({a: self(a) for a in alphabet
                             if a in self.domain().alphabet()})

    def _matrix_(self, R=None):
        r"""
        Return the incidence matrix of the morphism over the specified ring.

        EXAMPLES::

            sage: # needs sage.modules
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: tm = WordMorphism('a->ab,b->ba')
            sage: Mfibo = matrix(fibo); Mfibo      # indirect doctest
            [1 1]
            [1 0]
            sage: Mtm = matrix(tm); Mtm
            [1 1]
            [1 1]
            sage: Mtm * Mfibo == matrix(tm*fibo)   # indirect doctest
            True
            sage: Mfibo * Mtm == matrix(fibo*tm)   # indirect doctest
            True
            sage: Mfibo.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: p = Mfibo.charpoly(); p
            x^2 - x - 1
            sage: p.roots(ring=RR, multiplicities=False)
            [-0.618033988749895, 1.61803398874989]
        """
        if R is None:
            return self.incidence_matrix()
        else:
            return self.incidence_matrix().change_ring(R)

    def incidence_matrix(self):
        r"""
        Return the incidence matrix of the morphism. The order of the rows
        and column are given by the order defined on the alphabet of the
        domain and the codomain.

        The matrix returned is over the integers.  If a different ring is
        desired, use either the ``change_ring`` function or the ``matrix``
        function.

        EXAMPLES::

            sage: m = WordMorphism('a->abc,b->a,c->c')
            sage: m.incidence_matrix()                                                  # needs sage.modules
            [1 1 0]
            [1 0 0]
            [1 0 1]
            sage: m = WordMorphism('a->abc,b->a,c->c,d->abbccccabca,e->abc')
            sage: m.incidence_matrix()                                                  # needs sage.modules
            [1 1 0 3 1]
            [1 0 0 3 1]
            [1 0 1 5 1]
        """
        L = []
        domain_alphabet = self.domain().alphabet()
        codomain_alphabet = self.codomain().alphabet()
        for b in domain_alphabet:
            w = self._morph[b]
            ev_dict = w.evaluation_dict()
            L.append([ev_dict.get(a, 0) for a in codomain_alphabet])
        M = Matrix(IntegerRing(), L).transpose()
        return M

    def domain(self):
        r"""
        Return domain of ``self``.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').domain()
            Finite words over {'a', 'b'}
            sage: WordMorphism('b->ba,a->ab').domain()
            Finite words over {'a', 'b'}
            sage: WordMorphism('6->ab,y->5,0->asd').domain()
            Finite words over {'0', '6', 'y'}
        """
        return self._domain

    def codomain(self):
        r"""
        Return the codomain of ``self``.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').codomain()
            Finite words over {'a', 'b'}
            sage: WordMorphism('6->ab,y->5,0->asd').codomain()
            Finite words over {'5', 'a', 'b', 'd', 's'}
        """
        return self._codomain

    def is_endomorphism(self):
        r"""
        Return whether ``self`` is an endomorphism, that is if the
        domain coincide with the codomain.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').is_endomorphism()
            True
            sage: WordMorphism('6->ab,y->5,0->asd').is_endomorphism()
            False
            sage: WordMorphism('a->a,b->aa,c->aaa').is_endomorphism()
            False
            sage: Wabc = Words('abc')
            sage: m = WordMorphism('a->a,b->aa,c->aaa',codomain = Wabc)
            sage: m.is_endomorphism()
            True

        We check that :issue:`8674` is fixed::

            sage: P = WordPaths('abcd')                                                 # needs sage.modules
            sage: m = WordMorphism('a->adab,b->ab,c->cbcd,d->cd',                       # needs sage.modules
            ....:                  domain=P, codomain=P)
            sage: m.is_endomorphism()                                                   # needs sage.modules
            True
        """
        return self.codomain() == self.domain()

    def is_self_composable(self):
        r"""
        Return whether the codomain of ``self`` is contained in the domain.

        EXAMPLES::

            sage: f = WordMorphism('a->a,b->a')
            sage: f.is_endomorphism()
            False
            sage: f.is_self_composable()
            True
        """
        Adom = self.domain().alphabet()
        Acodom = self.codomain().alphabet()
        if Adom == Acodom:
            return True
        if Adom.cardinality() < Acodom.cardinality():
            return False
        if Adom.cardinality() == Infinity:
            raise NotImplementedError
        return all(a in Adom for a in Acodom)

    def image(self, letter):
        r"""
        Return the image of a letter.

        INPUT:

        - ``letter`` -- a letter in the domain alphabet

        OUTPUT: word

        .. NOTE::

            The letter is assumed to be in the domain alphabet
            (no check done). Hence, this method is faster
            than the ``__call__`` method suitable for words input.

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ac,c->a')
            sage: m.image('b')
            word: ac

        ::

            sage: s = WordMorphism({('a', 1):[('a', 1), ('a', 2)], ('a', 2):[('a', 1)]})
            sage: s.image(('a',1))
            word: ('a', 1),('a', 2)

        ::

            sage: s = WordMorphism({'b':[1,2], 'a':(2,3,4), 'z':[9,8,7]})
            sage: s.image('b')
            word: 12
            sage: s.image('a')
            word: 234
            sage: s.image('z')
            word: 987
        """
        return self._morph[letter]

    def images(self) -> list:
        r"""
        Return the list of all the images of the letters of the alphabet
        under ``self``.

        EXAMPLES::

            sage: sorted(WordMorphism('a->ab,b->a').images())
            [word: a, word: ab]
            sage: sorted(WordMorphism('6->ab,y->5,0->asd').images())
            [word: 5, word: ab, word: asd]
        """
        return list(self._morph.values())

    def reversal(self):
        r"""
        Return the reversal of ``self``.

        EXAMPLES::

            sage: WordMorphism('6->ab,y->5,0->asd').reversal()
            WordMorphism: 0->dsa, 6->ba, y->5
            sage: WordMorphism('a->ab,b->a').reversal()
            WordMorphism: a->ba, b->a
        """
        return WordMorphism({key: w.reversal()
                             for key, w in self._morph.items()},
                            codomain=self._codomain)

    def is_empty(self):
        r"""
        Return ``True`` if the cardinality of the domain is zero and
        ``False`` otherwise.

        EXAMPLES::

            sage: WordMorphism('').is_empty()
            True
            sage: WordMorphism('a->a').is_empty()
            False
        """
        return len(self._morph) == 0

    def is_erasing(self):
        r"""
        Return ``True`` if ``self`` is an erasing morphism, i.e. the image of a
        letter is the empty word.

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').is_erasing()
            False
            sage: WordMorphism('6->ab,y->5,0->asd').is_erasing()
            False
            sage: WordMorphism('6->ab,y->5,0->asd,7->').is_erasing()
            True
            sage: WordMorphism('').is_erasing()
            False
        """
        for image in self.images():
            if image.is_empty():
                return True
        return False

    def is_identity(self):
        r"""
        Return ``True`` if ``self`` is the identity morphism.

        EXAMPLES::

            sage: m = WordMorphism('a->a,b->b,c->c,d->e')
            sage: m.is_identity()
            False
            sage: WordMorphism('a->a,b->b,c->c').is_identity()
            True
            sage: WordMorphism('a->a,b->b,c->cb').is_identity()
            False
            sage: m = WordMorphism('a->b,b->c,c->a')
            sage: (m^2).is_identity()
            False
            sage: (m^3).is_identity()
            True
            sage: (m^4).is_identity()
            False
            sage: WordMorphism('').is_identity()
            True
            sage: WordMorphism({0:[0],1:[1]}).is_identity()
            True

        We check that :issue:`8618` is fixed::

            sage: t = WordMorphism({'a1':['a2'], 'a2':['a1']})
            sage: (t*t).is_identity()
            True
        """
        if not self.is_letter_permutation():
            return False
        for letter in self.domain().alphabet():
            img = self.image(letter)
            if img[0] != letter:
                return False
        return True

    def is_letter_permutation(self):
        r"""
        Return ``True`` if ``self`` is a letter permutation.

        EXAMPLES::

            sage: WordMorphism('a->b,b->a').is_letter_permutation()
            True
            sage: WordMorphism('a->b,b->ba').is_letter_permutation()
            False
            sage: WordMorphism({0:[1],1:[0]}).is_letter_permutation()
            True
            sage: WordMorphism('a->c,b->a,c->b').is_letter_permutation()
            True

        TESTS::

            sage: WordMorphism('').is_letter_permutation()
            True
            sage: WordMorphism({0:1,1:0,2:3}).is_letter_permutation()
            False
            sage: WordMorphism('a->c,b->c,c->b').is_letter_permutation()
            False
            sage: WordMorphism('a->').is_letter_permutation()
            False
        """
        if not self.is_endomorphism():
            return False
        for letter in self.domain().alphabet():
            img = self.image(letter)
            if img.length() != 1:
                return False
        return True

    def partition_of_domain_alphabet(self):
        r"""
        Return a partition of the domain alphabet.

        Let `\varphi:\Sigma^*\rightarrow\Sigma^*` be an involution. There
        exists a triple of sets `(A, B, C)` such that

         -  `A \cup B \cup C =\Sigma`;
         -  `A`, `B` and `C` are mutually disjoint and
         -  `\varphi(A)= B`, `\varphi(B)= A`, `\varphi(C)= C`.

        These sets are not unique.

        INPUT:

        - ``self`` -- an involution

        OUTPUT: a tuple of three sets

        EXAMPLES::

            sage: m = WordMorphism('a->b,b->a')
            sage: m.partition_of_domain_alphabet()  # random ordering
            ({'a'}, {'b'}, {})
            sage: m = WordMorphism('a->b,b->a,c->c')
            sage: m.partition_of_domain_alphabet()  # random ordering
            ({'a'}, {'b'}, {'c'})
            sage: m = WordMorphism('a->a,b->b,c->c')
            sage: m.partition_of_domain_alphabet()  # random ordering
            ({}, {}, {'a', 'c', 'b'})
            sage: m = WordMorphism('A->T,T->A,C->G,G->C')
            sage: m.partition_of_domain_alphabet()  # random ordering
            ({'A', 'C'}, {'T', 'G'}, {})
            sage: I = WordMorphism({0:oo,oo:0,1:-1,-1:1,2:-2,-2:2,3:-3,-3:3})
            sage: I.partition_of_domain_alphabet()  # random ordering
            ({0, -1, -3, -2}, {1, 2, 3, +Infinity}, {})

        TESTS::

            sage: m = WordMorphism('a->b,b->a,c->a')
            sage: m.partition_of_domain_alphabet()
            Traceback (most recent call last):
            ...
            TypeError: self is not an involution
        """
        if not self.is_involution():
            raise TypeError("self is not an involution")

        A = set()
        B = set()
        C = set()
        for a in self.domain().alphabet():
            if a == self(a)[0]:
                C.add(a)
            elif not (a in A or a in B):
                A.add(a)
                B.add(self(a)[0])

        return Set(A), Set(B), Set(C)

    def is_involution(self):
        r"""
        Return ``True`` if ``self`` is an involution, i.e. its square
        is the identity.

        INPUT:

        - ``self`` -- an endomorphism

        EXAMPLES::

            sage: WordMorphism('a->b,b->a').is_involution()
            True
            sage: WordMorphism('a->b,b->ba').is_involution()
            False
            sage: WordMorphism({0:[1],1:[0]}).is_involution()
            True
            sage: WordMorphism('a->,b->c').is_involution()
            False

        TESTS::

            sage: WordMorphism('').is_involution()
            True
            sage: WordMorphism({0:1,1:0,2:3}).is_involution()
            False
        """
        if not self.is_letter_permutation():
            return False
        return (self * self).is_identity()

    def pisot_eigenvector_right(self):
        r"""
        Return the right eigenvector of the incidence matrix associated
        to the largest eigenvalue (in absolute value).

        Unicity of the result is guaranteed when the multiplicity of the
        largest eigenvalue is one, for example when ``self`` is a Pisot
        irreductible substitution.

        A substitution is Pisot irreducible if the characteristic
        polynomial of its incidence matrix is irreducible over `\QQ` and
        has all roots, except one, of modulus strictly smaller than 1.

        INPUT:

        - ``self`` -- a Pisot irreducible substitution

        EXAMPLES::

            sage: m = WordMorphism('a->aaaabbc,b->aaabbc,c->aabc')
            sage: matrix(m)                                                             # needs sage.modules
            [4 3 2]
            [2 2 1]
            [1 1 1]
            sage: m.pisot_eigenvector_right()                                           # needs sage.modules sage.rings.number_field
            (1, 0.5436890126920763?, 0.2955977425220848?)
        """
        eig = self.incidence_matrix().eigenvectors_right()
        return max(eig, key=lambda x: abs(x[0]))[1][0]

    def pisot_eigenvector_left(self):
        r"""
        Return the left eigenvector of the incidence matrix associated
        to the largest eigenvalue (in absolute value).

        Unicity of the result is guaranteed when the multiplicity of the
        largest eigenvalue is one, for example when ``self`` is a Pisot
        irreductible substitution.

        A substitution is Pisot irreducible if the characteristic
        polynomial of its incidence matrix is irreducible over `\QQ` and
        has all roots, except one, of modulus strictly smaller than 1.

        INPUT:

        - ``self`` -- a Pisot irreducible substitution

        EXAMPLES::

            sage: m = WordMorphism('a->aaaabbc,b->aaabbc,c->aabc')
            sage: matrix(m)                                                             # needs sage.modules
            [4 3 2]
            [2 2 1]
            [1 1 1]
            sage: m.pisot_eigenvector_left()                                            # needs sage.modules sage.rings.number_field
            (1, 0.8392867552141611?, 0.5436890126920763?)
        """
        eig = self.incidence_matrix().eigenvectors_left()
        return max(eig, key=lambda x: abs(x[0]))[1][0]

    def _check_primitive(self):
        r"""
        Return ``True`` if all the letters of the domain appear in all the
        images of letters of the domain.

        INPUT:

        - ``self`` -- the codomain must be an instance of Words

        EXAMPLES::

            sage: m = WordMorphism('a->ab,b->ba')
            sage: m._check_primitive()
            True
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo._check_primitive()
            False
            sage: WordMorphism({2:[4,5,6],3:[4,1,8]})
            WordMorphism: 2->456, 3->418
            sage: WordMorphism({2:[4,5,6],3:[4,1,8]})._check_primitive()
            False
        """
        dom_alphabet = set(self.domain().alphabet())

        for image in self.images():
            if not dom_alphabet <= set(image):
                return False
        else:
            return True

    def is_primitive(self):
        r"""
        Return ``True`` if ``self`` is primitive.

        A morphism `\varphi` is *primitive* if there exists
        an positive integer `k` such that for all `\alpha\in\Sigma`,
        `\varphi^k(\alpha)` contains all the letters of `\Sigma`.

        INPUT:

        - ``self`` -- an endomorphism

        ALGORITHM:

            Exercices 8.7.8, p.281 in [1]:
            (c) Let `y(M)` be the least integer `e` such that `M^e` has all
            positive entries. Prove that, for all primitive matrices `M`,
            we have `y(M) \leq (d-1)^2 + 1`.
            (d) Prove that the bound `y(M)\leq (d-1)^2+1` is best possible.

        EXAMPLES::

            sage: tm = WordMorphism('a->ab,b->ba')
            sage: tm.is_primitive()                                                     # needs sage.modules
            True
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo.is_primitive()                                                   # needs sage.modules
            True
            sage: m = WordMorphism('a->bb,b->aa')
            sage: m.is_primitive()                                                      # needs sage.modules
            False
            sage: f = WordMorphism({0:[1],1:[0]})
            sage: f.is_primitive()                                                      # needs sage.modules
            False

        ::

            sage: s = WordMorphism('a->b,b->c,c->ab')
            sage: s.is_primitive()                                                      # needs sage.modules
            True
            sage: s = WordMorphism('a->b,b->c,c->d,d->e,e->f,f->g,g->h,h->ab')
            sage: s.is_primitive()                                                      # needs sage.modules
            True

        TESTS::

            sage: m = WordMorphism('a->bb,b->aac')
            sage: m.is_primitive()
            Traceback (most recent call last):
            ...
            TypeError: self (=a->bb, b->aac) is not an endomorphism
            sage: m = WordMorphism('a->,b->', codomain=Words('ab'))
            sage: m.is_primitive()                                                      # needs sage.modules
            False
            sage: m = WordMorphism('a->,b->')
            sage: m.is_primitive()
            Traceback (most recent call last):
            ...
            TypeError: self (=a->, b->) is not an endomorphism

        REFERENCES:

        - [1] Jean-Paul Allouche and Jeffrey Shallit, Automatic Sequences:
          Theory, Applications, Generalizations, Cambridge University Press,
          2003.
        """
        if not self.is_endomorphism():
            raise TypeError("self (=%s) is not an endomorphism" % self)
        return self.incidence_matrix().is_primitive()

    def is_prolongable(self, letter):
        r"""
        Return ``True`` if ``self`` is prolongable on ``letter``.

        A morphism `\varphi` is prolongable on a letter `a`
        if `a` is a prefix of `\varphi(a)`.

        INPUT:

        - ``self`` -- its codomain must be an instance of Words
        - ``letter`` -- a letter in the domain alphabet

        OUTPUT: boolean

        EXAMPLES::

            sage: WordMorphism('a->ab,b->a').is_prolongable(letter='a')
            True
            sage: WordMorphism('a->ab,b->a').is_prolongable(letter='b')
            False
            sage: WordMorphism('a->ba,b->ab').is_prolongable(letter='b')
            False
            sage: (WordMorphism('a->ba,b->ab')^2).is_prolongable(letter='b')
            True
            sage: WordMorphism('a->ba,b->').is_prolongable(letter='b')
            False
            sage: WordMorphism('a->bb,b->aac').is_prolongable(letter='a')
            False

        We check that :issue:`8595` is fixed::

            sage: s = WordMorphism({('a', 1) : [('a', 1), ('a', 2)], ('a', 2) : [('a', 1)]})
            sage: s.is_prolongable(('a',1))
            True

        TESTS::

            sage: WordMorphism('a->ab,b->b,c->ba').is_prolongable(letter='d')
            Traceback (most recent call last):
            ...
            TypeError: letter (=d) is not in the domain alphabet (={'a', 'b', 'c'})

        ::

            sage: n0, n1 = matrix(2,[1,1,1,0]), matrix(2,[2,1,1,0])                     # needs sage.modules
            sage: n = {'a':n0, 'b':n1}                                                  # needs sage.modules
            sage: WordMorphism(n).is_prolongable(letter='a')    # not implemented, needs sage.modules
            Traceback (most recent call last):
            ...
            TypeError: codomain of self must be an instance of Words
        """
        if letter not in self.domain().alphabet():
            raise TypeError("letter (={}) is not in the domain alphabet (={})".format(letter, self.domain().alphabet()))
        image = self.image(letter)
        return not image.is_empty() and letter == image[0]

    def is_uniform(self, k=None):
        r"""
        Return ``True`` if ``self`` is a `k`-uniform morphism.

        Let `k` be a positive integer. A morphism `\phi` is called `k`-uniform
        if for every letter `\alpha`, we have `|\phi(\alpha)| = k`. In other
        words, all images have length `k`. A morphism is called uniform if it
        is `k`-uniform for some positive integer `k`.

        INPUT:

        - ``k`` -- positive integer or ``None``. If set to a positive integer,
          then the function return ``True`` if ``self`` is `k`-uniform.
          If set to ``None``, then the function return ``True`` if ``self``
          is uniform.

        EXAMPLES::

            sage: phi = WordMorphism('a->ab,b->a')
            sage: phi.is_uniform()
            False
            sage: phi.is_uniform(k=1)
            False
            sage: tau = WordMorphism('a->ab,b->ba')
            sage: tau.is_uniform()
            True
            sage: tau.is_uniform(k=1)
            False
            sage: tau.is_uniform(k=2)
            True

        TESTS::

            sage: phi = WordMorphism('')
            sage: phi.is_uniform()
            True
        """
        if k is None:
            try:
                k = self.images()[0].length()
            except IndexError:
                return True
        return all(w.length() == k for w in self.images())

    def fixed_point(self, letter):
        r"""
        Return the fixed point of ``self`` beginning by the given ``letter``.

        A fixed point of morphism `\varphi` is a word `w` such that
        `\varphi(w) = w`.

        INPUT:

        - ``self`` -- an endomorphism (or more generally a self-composable
          morphism), must be prolongable on ``letter``

        - ``letter`` -- in the domain of ``self``, the first letter
          of the fixed point

        OUTPUT: ``word`` -- the fixed point of ``self`` beginning with ``letter``

        EXAMPLES::

            sage: W = FiniteWords('abc')

        1. Infinite fixed point::

            sage: WordMorphism('a->ab,b->ba').fixed_point(letter='a')
            word: abbabaabbaababbabaababbaabbabaabbaababba...
            sage: WordMorphism('a->ab,b->a').fixed_point(letter='a')
            word: abaababaabaababaababaabaababaabaababaaba...
            sage: WordMorphism('a->ab,b->b,c->ba', codomain=W).fixed_point(letter='a')
            word: abbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb...

        2. Infinite fixed point of an erasing morphism::

            sage: WordMorphism('a->ab,b->,c->ba', codomain=W).fixed_point(letter='a')
            word: ab

        3. Finite fixed point::

            sage: WordMorphism('a->ab,b->b,c->ba', codomain=W).fixed_point(letter='b')
            word: b
            sage: _.parent()
            Finite words over {'a', 'b', 'c'}

            sage: WordMorphism('a->ab,b->cc,c->', codomain=W).fixed_point(letter='a')
            word: abcc
            sage: _.parent()
            Finite words over {'a', 'b', 'c'}

            sage: m = WordMorphism('a->abc,b->,c->')
            sage: fp = m.fixed_point('a'); fp
            word: abc

            sage: m = WordMorphism('a->ba,b->')
            sage: m('ba')
            word: ba
            sage: m.fixed_point('a') #todo: not implemented
            word: ba

        5. Fixed point of a power of a morphism::

            sage: m = WordMorphism('a->ba,b->ab')
            sage: (m^2).fixed_point(letter='a')
            word: abbabaabbaababbabaababbaabbabaabbaababba...

        6. With a self-composable but not endomorphism

            sage: m = WordMorphism('a->cbc,b->bc,c->b')
            sage: m.is_endomorphism()
            False
            sage: m.fixed_point('b')
            word: bcbbcbcbbcbbcbcbbcbcbbcbbcbcbbcbbcbcbbcb...

        TESTS::

            sage: WordMorphism('a->ab,b->,c->ba', codomain=W).fixed_point(letter='b')
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on b
            sage: WordMorphism('a->ab,b->,c->ba', codomain=W).fixed_point(letter='c')
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on c
            sage: WordMorphism('a->ab,b->,c->ba', codomain=W).fixed_point(letter='d')
            Traceback (most recent call last):
            ...
            TypeError: letter (=d) is not in the domain alphabet (={'a', 'b', 'c'})
            sage: WordMorphism('a->aa,b->aac').fixed_point(letter='a')
            Traceback (most recent call last):
            ...
            TypeError: self (=a->aa, b->aac) is not self-composable
        """
        if not self.is_self_composable():
            raise TypeError("self (=%s) is not self-composable" % self)

        if not self.is_prolongable(letter=letter):
            raise TypeError("self must be prolongable on %s" % letter)

        parent = self.codomain()
        if self.is_growing(letter):
            from sage.combinat.words.word import InfiniteWord_morphic
            return InfiniteWord_morphic(parent.shift(), self, letter,
                                        coding=None, length=Infinity)
        else:
            from sage.combinat.words.word import FiniteWord_morphic
            w = FiniteWord_morphic(parent, self, letter,
                                   coding=None, length='finite')
            # since FiniteWord_morphic uses the method __getitem__
            # from FiniteWord_callable, the length must be precomputed
            # for __getitem__ to work properly
            w.length()
            return w

    def fixed_points(self):
        r"""
        Return the list of all fixed points of ``self``.

        EXAMPLES::

            sage: f = WordMorphism('a->ab,b->ba')
            sage: for w in f.fixed_points(): print(w)
            abbabaabbaababbabaababbaabbabaabbaababba...
            baababbaabbabaababbabaabbaababbaabbabaab...

            sage: f = WordMorphism('a->ab,b->c,c->a')
            sage: for w in f.fixed_points(): print(w)
            abcaababcabcaabcaababcaababcabcaababcabc...

            sage: f = WordMorphism('a->ab,b->cab,c->bcc')
            sage: for w in f.fixed_points(): print(w)
            abcabbccabcabcabbccbccabcabbccabcabbccab...

        This shows that issue :issue:`13668` has been resolved::

            sage: d = {1:[1,2],2:[2,3],3:[4],4:[5],5:[6],6:[7],7:[8],8:[9],9:[10],10:[1]}
            sage: s = WordMorphism(d)
            sage: s7 = s^7
            sage: s7.fixed_points()
            [word: 12232342..., word: 2,3,4,5,6,7,8...]
            sage: s7r = s7.reversal()
            sage: s7r.periodic_point(2)
            word: 2,1,1,10,9,8,7,6,5,4,3,2,1,10,9,8,7,6,5,4,3,2,10,9,8,7,6,5,4,3,2,9,8,7,6,5,4,3,2,8,...

        This shows that issue :issue:`13668` has been resolved::

            sage: s = "1->321331332133133,2->133321331332133133,3->2133133133321331332133133"
            sage: s = WordMorphism(s)
            sage: (s^2).fixed_points()
            []
        """
        return [self.fixed_point(letter=letter)
                for letter in self.domain().alphabet()
                if self.is_prolongable(letter=letter)]

    def periodic_point(self, letter):
        r"""
        Return the periodic point of ``self`` that starts with ``letter``.

        EXAMPLES::

            sage: f = WordMorphism('a->bab,b->ab')
            sage: f.periodic_point('a')
            word: abbababbababbabababbababbabababbababbaba...
            sage: f.fixed_point('a')
            Traceback (most recent call last):
            ...
            TypeError: self must be prolongable on a

        Make sure that :issue:`31759` is fixed::

            sage: WordMorphism('a->b,b->a').periodic_point('a')
            word: a
        """
        if not self.is_growing(letter):
            w = self.domain()(letter)
            prev = set()
            while w not in prev:
                prev.add(w)
                w = self(w)
            return w

        elif self.is_erasing():
            raise NotImplementedError("self should be non erasing")

        else:
            cycle = [letter]
            a = self(letter)[0]
            while a not in cycle:
                cycle.append(a)
                a = self(a)[0]
            if a != letter:
                raise ValueError("there is no periodic point starting with letter (=%s)" % letter)

            P = PeriodicPointIterator(self, cycle)
            return self.codomain().shift()(P._cache[0])

    def periodic_points(self):
        r"""
        Return the periodic points of ``f`` as a list of tuples where each tuple is
        a periodic orbit of ``f``.

        EXAMPLES::

            sage: f = WordMorphism('a->aba,b->baa')
            sage: for p in f.periodic_points():
            ....:     print("{} , {}".format(len(p), p[0]))
            1 , ababaaababaaabaabaababaaababaaabaabaabab...
            1 , baaabaabaababaaabaababaaabaababaaababaaa...

            sage: f = WordMorphism('a->bab,b->aa')
            sage: for p in f.periodic_points():
            ....:     print("{} , {}".format(len(p), p[0]))
            2 , aababaaaababaababbabaababaababbabaababaa...
            sage: f.fixed_points()
            []

        This shows that issue :issue:`13668` has been resolved::

            sage: d = {1:[1,2],2:[2,3],3:[4],4:[5],5:[6],6:[7],7:[8],8:[9],9:[10],10:[1]}
            sage: s = WordMorphism(d)
            sage: s7 = s^7
            sage: s7r = s7.reversal()
            sage: for p in s7r.periodic_points(): p
            [word: 1,10,9,8,7,6,5,4,3,2,10,9,8,7,6,5,4,3,2,...,
             word: 8765432765432654325432432322176543265432...,
             word: 5,4,3,2,4,3,2,3,2,2,1,4,3,2,3,2,2,1,3,2,...,
             word: 2,1,1,10,9,8,7,6,5,4,3,2,1,10,9,8,7,6,5,...,
             word: 9876543287654327654326543254324323221876...,
             word: 6543254324323221543243232214323221322121...,
             word: 3,2,2,1,2,1,1,10,9,8,7,6,5,4,3,2,2,1,1,1...,
             word: 10,9,8,7,6,5,4,3,2,9,8,7,6,5,4,3,2,8,7,6...,
             word: 7654326543254324323221654325432432322154...,
             word: 4,3,2,3,2,2,1,3,2,2,1,2,1,1,10,9,8,7,6,5...]

        Make sure that :issue:`31454` is fixed::

            sage: WordMorphism('a->a,b->bb').periodic_points()
            [[word: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb...]]
        """
        if not self.is_endomorphism():
            raise ValueError("f should be an endomorphism")

        if self.is_erasing():
            raise NotImplementedError("f should be non erasing")

        A = self.domain().alphabet()
        d = {letter: self(letter)[0] for letter in A}
        G = set(self.growing_letters())

        res = []
        parent = self.codomain().shift()
        for cycle in get_cycles(CallableDict(d), A):
            if cycle[0] in G:
                P = PeriodicPointIterator(self, cycle)
                res.append([parent(P._cache[i]) for i in range(len(cycle))])

        return res

    def _language_naive(self, n, u):
        r"""
        Return all words of length less than ``n`` by naive substitution.

        The language of the substitution is the DOL language which consist
        of factors of `s^n(u)`.

        This method assumes this substitution is non-erasing.

        INPUT:

        - ``n`` -- nonnegative integer; length of the words in the language

        - ``u`` -- a word used as a seed

        OUTPUT: a Python set

        TESTS::

            sage: s = WordMorphism({0: [0,1], 1:[0]})
            sage: W = s.domain()
            sage: sorted(s._language_naive(3, W([0])))
            [word: 0, word: 00, word: 01, word: 1, word: 10]
            sage: sorted(s._language_naive(3, W([1])))
            [word: 0, word: 00, word: 01, word: 1, word: 10]

            sage: s._language_naive(3, W())
            set()
            sage: W([1, 1]) in s._language_naive(3, W([1, 1]))
            True
        """
        L = set()
        todo = []
        for i in range(len(u)):
            for j in range(i + 1, min(len(u) + 1, i + n)):
                f = u[i:j]
                if f not in L:
                    todo.append(f)
                    L.add(f)
        while todo:
            u = todo.pop()
            v = self(u)
            if u.length() == 1:
                for i in range(len(v)):
                    for j in range(i + 1, min(len(v) + 1, i + n)):
                        f = v[i:j]
                        if f not in L:
                            todo.append(f)
                            L.add(f)
            else:
                l = self._morph[u[0]].length()
                r = self._morph[u[-1]].length()
                m = v.length() - l - r
                x = n - 1 - m
                for i in range(l - min(x - 1, l), l):
                    for j in range(l + m + 1, l + m + 1 + min(x - l + i, r)):
                        f = v[i:j]
                        if f not in L:
                            todo.append(f)
                            L.add(f)
        return L

    def language(self, n, u=None):
        r"""
        Return the words of length ``n`` in the language generated by this substitution.

        Given a non-erasing substitution `s` and a word `u` the DOL-language
        generated by `s` and `u` is the union of the factors of `s^n(u)` where
        `n` is a nonnegative integer.

        INPUT:

        - ``n`` -- nonnegative integer; length of the words in the language

        - ``u`` -- a word or ``None`` (default: ``None``); if set to
          ``None`` some letter of the alphabet is used

        OUTPUT: a Python set

        EXAMPLES:

        The fibonacci morphism::

            sage: s = WordMorphism({0: [0,1], 1: [0]})
            sage: sorted(s.language(3))                                                 # needs sage.modules
            [word: 001, word: 010, word: 100, word: 101]
            sage: len(s.language(1000))                                                 # needs sage.modules
            1001
            sage: all(len(s.language(n)) == n+1 for n in range(100))                    # needs sage.modules
            True

        A growing but non-primitive example. The DOL-languages generated
        by 0 and 2 are different::

            sage: s = WordMorphism({0: [0,1], 1: [0], 2: [2,0,2]})

            sage: u = s.fixed_point(0)
            sage: A0 = u[:200].factor_set(5)                                            # needs sage.modules
            sage: B0 = s.language(5, [0])                                               # needs sage.modules
            sage: set(A0) == B0                                                         # needs sage.modules
            True

            sage: v = s.fixed_point(2)
            sage: A2 = v[:200].factor_set(5)                                            # needs sage.modules
            sage: B2 = s.language(5, [2])                                               # needs sage.modules
            sage: set(A2) == B2                                                         # needs sage.modules
            True

            sage: len(A0), len(A2)                                                      # needs sage.modules
            (6, 20)

        The Chacon transformation (non-primitive)::

            sage: s = WordMorphism({0: [0,0,1,0], 1:[1]})
            sage: sorted(s.language(10))                                                # needs sage.modules
            [word: 0001000101,
             word: 0001010010,
             ...
             word: 1010010001,
             word: 1010010100]
        """
        W = self.domain()
        if self.codomain() != W:
            raise ValueError('substitution not an endomorphism')

        if n == 0:
            return [W()]

        A = W.alphabet()
        if u is None:
            u = W([A.an_element()])
        else:
            u = W(u)

        if n <= 2 or not self.is_growing():
            return [w for w in self._language_naive(n + 1, u) if len(w) == n]

        # compute the right power
        M = m = self.incidence_matrix().transpose()
        p = 1
        d = m.nrows()
        while any(sum(M.row(j)) < n for j in range(d)):
            M *= m
            p += 1
        s = self**p
        im = {a: s.image(a) for a in A}

        # build factors by considering concatenations of images
        # of two letter words
        L2 = (w for w in self._language_naive(3, u) if len(w) == 2)
        L = set()
        for u in L2:
            v = im[u[0]] + im[u[1]]
            for k in range(len(v) - n + 1):
                L.add(v[k:k + n])
        return L

    def conjugate(self, pos):
        r"""
        Return the morphism where the image of the letter by ``self``
        is conjugated of parameter ``pos``.

        INPUT:

        - ``pos`` -- integer

        EXAMPLES::

            sage: m = WordMorphism('a->abcde')
            sage: m.conjugate(0) == m
            True
            sage: m.conjugate(1)
            WordMorphism: a->bcdea
            sage: m.conjugate(3)
            WordMorphism: a->deabc
            sage: WordMorphism('').conjugate(4)
            WordMorphism:
            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: m.conjugate(2)
            WordMorphism: a->cdeab, b->zxy
        """
        return WordMorphism({key: w.conjugate(pos)
                             for (key, w) in self._morph.items()})

    def has_left_conjugate(self):
        r"""
        Return ``True`` if all the non empty images of ``self`` begins with
        the same letter.

        EXAMPLES::

            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: m.has_left_conjugate()
            False
            sage: WordMorphism('b->xyz').has_left_conjugate()
            True
            sage: WordMorphism('').has_left_conjugate()
            True
            sage: WordMorphism('a->,b->xyz').has_left_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb').has_left_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb,c->').has_left_conjugate()
            True
        """
        I = (w for w in self.images() if not FiniteWord_class.is_empty(w))

        try:
            letter = next(I)[0]
        except StopIteration:
            return True

        # Compare the first letter of all the non empty images
        return all(image[0] == letter for image in I)

    def has_right_conjugate(self):
        r"""
        Return ``True`` if all the non empty images of ``self`` ends with the
        same letter.

        EXAMPLES::

            sage: m = WordMorphism('a->abcde,b->xyz')
            sage: m.has_right_conjugate()
            False
            sage: WordMorphism('b->xyz').has_right_conjugate()
            True
            sage: WordMorphism('').has_right_conjugate()
            True
            sage: WordMorphism('a->,b->xyz').has_right_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb').has_right_conjugate()
            True
            sage: WordMorphism('a->abbab,b->abb,c->').has_right_conjugate()
            True
        """
        return self.reversal().has_left_conjugate()

    def list_of_conjugates(self):
        r"""
        Return the list of all the conjugate morphisms of ``self``.

        DEFINITION:

        Recall from Lothaire [1] (Section 2.3.4)
        that `\varphi` is *right conjugate* of `\varphi'`,
        noted `\varphi\triangleleft\varphi'`, if there exists
        `u \in \Sigma^*` such that

        .. MATH::

            \varphi(\alpha)u = u\varphi'(\alpha),

        for all `\alpha \in \Sigma`, or equivalently that
        `\varphi(x)u = u\varphi'(x)`, for all words `x \in \Sigma^*`.
        Clearly, this relation is not
        symmetric so that we say that two morphisms `\varphi` and
        `\varphi'` are *conjugate*, noted
        `\varphi\bowtie\varphi'`, if
        `\varphi\triangleleft\varphi'` or
        `\varphi'\triangleleft\varphi`. It is easy to see that
        conjugacy of morphisms is an equivalence relation.

        REFERENCES:

        - [1] M. Lothaire, Algebraic Combinatorics on words, Cambridge
          University Press, 2002.

        EXAMPLES::

            sage: m = WordMorphism('a->abbab,b->abb')
            sage: m.list_of_conjugates()
            [WordMorphism: a->babba, b->bab,
            WordMorphism: a->abbab, b->abb,
            WordMorphism: a->bbaba, b->bba,
            WordMorphism: a->babab, b->bab,
            WordMorphism: a->ababb, b->abb,
            WordMorphism: a->babba, b->bba,
            WordMorphism: a->abbab, b->bab]
            sage: m = WordMorphism('a->aaa,b->aa')
            sage: m.list_of_conjugates()
            [WordMorphism: a->aaa, b->aa]
            sage: WordMorphism('').list_of_conjugates()
            [WordMorphism: ]
            sage: m = WordMorphism('a->aba,b->aba')
            sage: m.list_of_conjugates()
            [WordMorphism: a->baa, b->baa,
            WordMorphism: a->aab, b->aab,
            WordMorphism: a->aba, b->aba]
            sage: m = WordMorphism('a->abb,b->abbab,c->')
            sage: m.list_of_conjugates()
            [WordMorphism: a->bab, b->babba, c->,
            WordMorphism: a->abb, b->abbab, c->,
            WordMorphism: a->bba, b->bbaba, c->,
            WordMorphism: a->bab, b->babab, c->,
            WordMorphism: a->abb, b->ababb, c->,
            WordMorphism: a->bba, b->babba, c->,
            WordMorphism: a->bab, b->abbab, c->]
        """
        if self.is_empty():
            return [self]

        # Build the list c of conjugate morphisms
        c = []
        m = self
        c.append(m)
        while m.has_left_conjugate():
            m = m.conjugate(1)
            if m == self:
                break
            c.append(m)
        m = self
        while m.has_right_conjugate():
            m = m.conjugate(-1)
            if m == self:
                break
            c.insert(0, m)

        # Build the list d of distinct morphisms
        d = []
        for m in c:
            if m not in d:
                d.append(m)
        return d

    def is_in_classP(self, f=None):
        r"""
        Return ``True`` if ``self`` is in class `P` (or `f`-`P`).

        DEFINITION : Let `A` be an alphabet. We say that a
        primitive substitution `S` is in the *class P* if there
        exists a palindrome `p` and for each `b\in A` a
        palindrome `q_b` such that `S(b)=pq_b` for all
        `b\in A`. [1]

        Let `f` be an involution on `A`. "We say that a morphism
        `\varphi` is in class `f`-`P` if there exists an
        `f`-palindrome `p` and for each `\alpha \in A`
        there exists an `f`-palindrome `q_\alpha` such
        that `\varphi(\alpha)=pq_\alpha`. [2]

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. WordMorphism)

        REFERENCES:

        - [1] Hof, A., O. Knill et B. Simon, Singular continuous
          spectrum for palindromic Schrödinger operators,
          Commun. Math. Phys.  174 (1995) 149-159.

        - [2] Labbe, Sebastien. Proprietes combinatoires des
          `f`-palindromes, Memoire de maitrise en Mathematiques,
          Montreal, UQAM, 2008, 109 pages.

        EXAMPLES::

            sage: WordMorphism('a->bbaba,b->bba').is_in_classP()
            True
            sage: tm = WordMorphism('a->ab,b->ba')
            sage: tm.is_in_classP()
            False
            sage: f = WordMorphism('a->b,b->a')
            sage: tm.is_in_classP(f=f)
            True
            sage: (tm^2).is_in_classP()
            True
            sage: (tm^2).is_in_classP(f=f)
            False
            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo.is_in_classP()
            True
            sage: fibo.is_in_classP(f=f)
            False
            sage: (fibo^2).is_in_classP()
            False
            sage: f = WordMorphism('a->b,b->a,c->c')
            sage: WordMorphism('a->acbcc,b->acbab,c->acbba').is_in_classP(f)
            True
        """
        if self.is_empty():
            return True

        # Compute the longest common prefix of all the images of letters
        images = self.images()
        lcp = images[0]
        for image in images:
            lcp = lcp.longest_common_prefix(image)

        # Find a common palindrome prefix
        for i in range(lcp.length() + 1):
            if lcp[:i].is_palindrome(f=f):

                # If all the suffixes are palindromes,
                for image in images:
                    if not image[i:].is_palindrome(f=f):
                        break
                else:
                    return True

        return False

    def has_conjugate_in_classP(self, f=None):
        r"""
        Return ``True`` if ``self`` has a conjugate in class `f`-`P`.

        DEFINITION : Let `A` be an alphabet. We say that a
        primitive substitution `S` is in the *class P* if there
        exists a palindrome `p` and for each `b\in A` a
        palindrome `q_b` such that `S(b)=pq_b` for all
        `b\in A`. [1]

        Let `f` be an involution on `A`. We say that a morphism
        `\varphi` is in class `f`-`P` if there exists an
        `f`-palindrome `p` and for each `\alpha \in A`
        there exists an `f`-palindrome `q_\alpha` such
        that `\varphi(\alpha)=pq_\alpha`. [2]

        INPUT:

        - ``f`` -- involution (default: ``None``) on the alphabet of ``self``;
          it must be callable on letters as well as words (e.g. WordMorphism)

        REFERENCES:

        - [1] Hof, A., O. Knill et B. Simon, Singular continuous
          spectrum for palindromic Schrödinger operators,
          Commun. Math. Phys.  174 (1995) 149-159.

        - [2] Labbe, Sebastien. Proprietes combinatoires des
          `f`-palindromes, Memoire de maitrise en Mathematiques,
          Montreal, UQAM, 2008, 109 pages.

        EXAMPLES::

            sage: fibo = WordMorphism('a->ab,b->a')
            sage: fibo.has_conjugate_in_classP()
            True
            sage: (fibo^2).is_in_classP()
            False
            sage: (fibo^2).has_conjugate_in_classP()
            True
        """
        for k in self.list_of_conjugates():
            if k.is_in_classP(f=f):
                return True
        return False

    def dual_map(self, k=1):
        r"""
        Return the dual map `E_k^*` of ``self`` (see [1]).

        .. NOTE::

            It is actually implemented only for `k=1`.

        INPUT:

        - ``self`` -- unimodular endomorphism defined on integers
          ``1, 2, \ldots, d``
        - ``k`` -- integer (default: 1)

        OUTPUT: an instance of E1Star - the dual map

        EXAMPLES::

            sage: sigma = WordMorphism({1: [2], 2: [3], 3: [1,2]})
            sage: sigma.dual_map()                                                      # needs sage.modules
            E_1^*(1->2, 2->3, 3->12)

        ::

            sage: sigma.dual_map(k=2)
            Traceback (most recent call last):
            ...
            NotImplementedError: the dual map E_k^* is implemented only for k = 1 (not 2)

        REFERENCES:

        - [1] Sano, Y., Arnoux, P. and Ito, S., Higher dimensional
          extensions of substitutions and their dual maps, Journal
          d'Analyse Mathématique 83 (2001), 183-206.
        """
        if k == 1:
            from sage.combinat.e_one_star import E1Star
            return E1Star(self)

        raise NotImplementedError("the dual map E_k^* is implemented only "
                                  "for k = 1 (not %s)" % k)

    @cached_method
    def rauzy_fractal_projection(self, eig=None, prec=53):
        r"""
        Return a dictionary giving the projection of the canonical basis.

        See the method :meth:`rauzy_fractal_plot` for more details about the projection.

        INPUT:

        - ``eig`` -- a real element of ``QQbar`` of degree >= 2 (default: ``None``).
          The eigenvalue used for the projection.
          It must be an eigenvalue of ``self.incidence_matrix()``.
          The one used by default is the maximal eigenvalue of
          ``self.incidence_matrix()`` (usually a Pisot number),
          but for substitutions with more than 3 letters
          other interesting choices are sometimes possible.

        - ``prec`` -- integer (default: 53);
          the number of bits used in the floating point representations
          of the coordinates

        OUTPUT:

            dictionary, letter -> vector, giving the projection

        EXAMPLES:

        The projection for the Rauzy fractal of the Tribonacci substitution
        is::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: s.rauzy_fractal_projection()                                          # needs sage.modules
            {'1': (1.00000000000000, 0.000000000000000),
             '2': (-1.41964337760708, -0.606290729207199),
             '3': (-0.771844506346038, 1.11514250803994)}

        TESTS::

            sage: t = WordMorphism('1->12,2->3,3->45,4->5,5->6,6->7,7->8,8->1')
            sage: E = t.incidence_matrix().eigenvalues()                                # needs sage.modules
            sage: x = [x for x in E if -0.8 < x < -0.7][0]                              # needs sage.modules
            sage: t.rauzy_fractal_projection(prec=10)                                   # needs sage.modules
            {'1': (1.0, 0.00),
             '2': (-1.7, -0.56),
             '3': (0.79, 1.3),
             '4': (1.9, -0.74),
             '5': (-1.7, -0.56),
             '6': (0.79, 1.3),
             '7': (0.21, -1.3),
             '8': (-0.88, 0.74)}
            sage: t.rauzy_fractal_projection(eig=x, prec=10)                            # needs sage.modules
            {'1': (1.0, 0.00),
             '2': (-0.12, -0.74),
             '3': (-0.66, -0.56),
             '4': (-0.46, -0.18),
             '5': (-0.54, 0.18),
             '6': (-0.34, 0.56),
             '7': (0.12, 0.74),
             '8': (0.66, 0.56)}

        AUTHOR:

            Timo Jolivet (2012-06-16)
        """
        alphabet = self.domain().alphabet()
        size_alphabet = len(alphabet)

        # Eigenvalues
        if eig is None:
            beta = max(self.incidence_matrix().eigenvalues(), key=abs)
        else:
            beta = eig

        # Test is deg(beta) >= 2
        if beta.degree() < 2:
            raise ValueError("the algebraic degree of ``eig`` must be at least two")

        # Algebraic conjugates of beta
        from sage.rings.qqbar import QQbar
        beta_conjugates = beta.minpoly().roots(QQbar, multiplicities=False)
        if not beta.imag():
            beta_conjugates.remove(beta)
        for x in beta_conjugates:
            if x.imag():
                beta_conjugates.remove(x.conjugate())

        # Left eigenvector vb in the number field Q(beta)
        from sage.rings.number_field.number_field import NumberField
        K = NumberField(beta.minpoly(), 'b')
        vb = (self.incidence_matrix() - K.gen()).kernel().basis()[0]

        # Projections of canonical base vectors from R^size_alphabet to C, using vb
        from sage.modules.free_module import VectorSpace
        canonical_basis = VectorSpace(K, size_alphabet).basis()
        canonical_basis_proj = {}

        from sage.rings.real_mpfr import RealField
        RealField_prec = RealField(prec)
        for a, x in zip(alphabet, canonical_basis):
            v = []
            for y in beta_conjugates:
                # if y has nonzero imaginary part
                if y.imag():
                    z = (vb * x).lift()(y)
                    z1, z2 = z.real(), z.imag()
                    v += [RealField_prec(z1), RealField_prec(z2)]
                # if y is real
                else:
                    z = (vb * x).lift()(y)
                    v += [RealField_prec(z)]
            canonical_basis_proj[a] = vector(v)

        return canonical_basis_proj

    def rauzy_fractal_points(self, n=None, exchange=False, eig=None, translate=None, prec=53):
        r"""
        Return a dictionary of list of points associated with the pieces
        of the Rauzy fractal of ``self``.

        INPUT:

        See the method :meth:`rauzy_fractal_plot` for a description
        of the options and more examples.

        OUTPUT: dictionary of list of points

        EXAMPLES:

        The Rauzy fractal of the Tribonacci substitution and the number of
        points in the piece of the fractal associated with ``'1'``, ``'2'``
        and ``'3'`` are respectively::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: D = s.rauzy_fractal_points(n=100)                                     # needs sage.modules
            sage: len(D['1'])                                                           # needs sage.modules
            54
            sage: len(D['2'])                                                           # needs sage.modules
            30
            sage: len(D['3'])                                                           # needs sage.modules
            16

        TESTS::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: D = s.rauzy_fractal_points(n=100, exchange=True,                      # needs sage.modules
            ....:                            translate=[(3,1,-2), (5,-33,8)], prec=40)
            sage: len(D['1'])                                                           # needs sage.modules
            108

        AUTHOR:

            Timo Jolivet (2012-06-16)
        """
        alphabet = self.domain().alphabet()
        canonical_basis_proj = self.rauzy_fractal_projection(eig=eig, prec=prec)

        # if exchange, set the projection to its opposite
        if exchange:
            for a in canonical_basis_proj:
                canonical_basis_proj[a] = - canonical_basis_proj[a]

        # Compute a fixed point u
        if exchange:
            u = iter(self.reversal().periodic_points()[0][0])
        else:
            u = iter(self.periodic_points()[0][0])

        # Manage various options in function of dimension
        if n is None:
            dim_fractal = len(canonical_basis_proj[alphabet[0]])
            if dim_fractal == 1:
                n = 1000
            elif dim_fractal == 2:
                n = 50000
            elif dim_fractal == 3:
                n = 5000
            else:
                n = 50000

        # Compute orbit points to plot
        S = 0
        orbit_points = {a: [] for a in alphabet}
        for _ in range(n):
            a = next(u)
            S += canonical_basis_proj[a]
            orbit_points[a].append(S)

        # Manage translated copies
        from sage.rings.real_mpfr import RealField
        RealField_prec = RealField(prec)
        if translate is not None:

            if isinstance(translate, dict):
                for a in translate:
                    translate[a] = [vector(RealField_prec, v) for v in translate[a]]

            else:
                translate = [vector(RealField_prec, v) for v in translate]

            for a in alphabet:
                translated_copies = {i: [] for i in alphabet}

                if isinstance(translate, list):
                    to_treat = translate

                elif isinstance(translate, dict):
                    try:
                        to_treat = translate[a]
                    except KeyError:
                        to_treat = []

                for x in to_treat:
                    v = 0
                    for i, z in zip(alphabet, x):
                        v += z * canonical_basis_proj[i]
                    translated_copies[a] += [vector(v) + w for w in orbit_points[a]]

                orbit_points[a] = translated_copies[a]

        return orbit_points

    def rauzy_fractal_plot(self, n=None, exchange=False, eig=None,
                           translate=None, prec=53,
                           colormap='hsv', opacity=None, plot_origin=None,
                           plot_basis=False, point_size=None):
        r"""
        Return a plot of the Rauzy fractal associated with a substitution.

        The substitution does not have to be irreducible.
        The usual definition of a Rauzy fractal requires that
        its dominant eigenvalue is a Pisot number but the present method
        doesn't require this, allowing to plot some interesting pictures
        in the non-Pisot case (see the examples below).

        For more details about the definition of the fractal and the
        projection which is used, see Section 3.1 of [1].

        Plots with less than 100,000 points take a few seconds,
        and several millions of points can be plotted in reasonable time.

        Other ways to draw Rauzy fractals (and more generally projections of paths)
        can be found in :meth:`sage.combinat.words.paths.FiniteWordPath_all.plot_projection`
        or in :meth:`sage.combinat.e_one_star`.

        OUTPUT: a Graphics object

        INPUT:

        - ``n`` -- integer (default: ``None``)
          The number of points used to plot the fractal.
          Default values: ``1000`` for a 1D fractal,
          ``50000`` for a 2D fractal, ``10000`` for a 3D fractal.

        - ``exchange`` -- boolean (default: ``False``); plot the Rauzy fractal
          with domain exchange

        - ``eig`` -- a real element of ``QQbar`` of degree >= 2 (default: ``None``);
          the eigenvalue used to plot the fractal.
          It must be an eigenvalue of ``self.incidence_matrix()``.
          The one used by default the maximal eigenvalue of
          ``self.incidence_matrix()`` (usually a Pisot number),
          but for substitutions with more than 3 letters
          other interesting choices are sometimes possible.

        - ``translate`` -- list of vectors of ``RR^size_alphabet``,
          or a dictionary from the alphabet to lists of vectors (default: ``None``).
          Plot translated copies of the fractal.
          This option allows to plot tilings easily.
          The projection used for these vectors is the same as
          the projection used for the canonical basis to plot the fractal.
          If the input is a list, all the pieces will be translated and plotted.
          If the input is a dictionary, each piece will be translated and plotted
          accordingly to the vectors associated with each letter in the dictionary.
          Note: by default, the Rauzy fractal placed at the origin
          is not plotted with the ``translate`` option;
          the vector ``(0,0,...,0)`` has to be added manually.

        - ``prec`` -- integer (default: 53);
          the number of bits used in the floating point representations
          of the points of the fractal

        - ``colormap`` -- color map or dictionary (default: ``'hsv'``).
          It can be one of the following:

           - ``string`` -- a coloring map. For available coloring map names type:
             ``sorted(colormaps)``

           - ``dict`` -- dictionary of the alphabet mapped to colors

        - ``opacity`` -- dictionary from the alphabet to the real interval
          [0,1] (default: ``None``); if none is specified, all letters are
          plotted with opacity ``1``

        - ``plot_origin`` -- a couple ``(k,c)`` (default: ``None``);
          if specified, mark the origin by a point of size ``k`` and color ``c``

        - ``plot_basis`` -- boolean (default: ``False``); plot the projection
          of the canonical basis with the fractal

        - ``point_size`` -- float (default: ``None``); the size of the points
          used to plot the fractal

        EXAMPLES:

        #. The Rauzy fractal of the Tribonacci substitution::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: s.rauzy_fractal_plot()        # long time                             # needs sage.plot
            Graphics object consisting of 3 graphics primitives

        #. The "Hokkaido" fractal. We tweak the plot using the plotting options
           to get a nice reusable picture, in which we mark the origin by a black dot::

            sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
            sage: G = s.rauzy_fractal_plot(n=100000, point_size=3,          # not tested
            ....:                          plot_origin=(50,"black"))
            sage: G.show(figsize=10, axes=false) # not tested

        #. Another "Hokkaido" fractal and its domain exchange::

            sage: s = WordMorphism({1:[2], 2:[4,3], 3:[4], 4:[5,3], 5:[6], 6:[1]})
            sage: s.rauzy_fractal_plot()                                    # not tested (> 1 second)
            sage: s.rauzy_fractal_plot(exchange=True)                       # not tested (> 1 second)

        #. A three-dimensional Rauzy fractal::

            sage: s = WordMorphism('1->12,2->13,3->14,4->1')
            sage: s.rauzy_fractal_plot()                                    # not tested (> 1 second)

        #. A one-dimensional Rauzy fractal (very scattered)::

            sage: s = WordMorphism('1->2122,2->1')
            sage: s.rauzy_fractal_plot().show(figsize=20)                   # not tested (> 1 second)

        #. A high resolution plot of a complicated fractal::

            sage: s = WordMorphism('1->23,2->123,3->1122233')
            sage: G = s.rauzy_fractal_plot(n=300000)                        # not tested (> 1 second)
            sage: G.show(axes=false, figsize=20)                            # not tested (> 1 second)

        #. A nice colorful animation of a domain exchange::

            sage: s = WordMorphism('1->21,2->3,3->4,4->25,5->6,6->7,7->1')
            sage: L = [s.rauzy_fractal_plot(),                              # not tested (> 1 second)
            ....:      s.rauzy_fractal_plot(exchange=True)]
            sage: animate(L, axes=false).show(delay=100)                    # not tested (> 1 second)

        #. Plotting with only one color::

            sage: s = WordMorphism('1->12,2->31,3->1')
            sage: cm = {'1':'black', '2':'black', '3':'black'}
            sage: s.rauzy_fractal_plot(colormap=cm)                         # not tested (> 1 second)

        #. Different fractals can be obtained by choosing another (non-Pisot) eigenvalue::

            sage: s = WordMorphism('1->12,2->3,3->45,4->5,5->6,6->7,7->8,8->1')
            sage: E = s.incidence_matrix().eigenvalues()                                # needs sage.modules
            sage: x = [x for x in E if -0.8 < x < -0.7][0]                              # needs sage.modules
            sage: s.rauzy_fractal_plot()                                    # not tested (> 1 second)
            sage: s.rauzy_fractal_plot(eig=x)                               # not tested (> 1 second)

        #. A Pisot reducible substitution with seemingly overlapping tiles::

            sage: s = WordMorphism({1:[1,2], 2:[2,3], 3:[4], 4:[5], 5:[6],
            ....:                   6:[7], 7:[8], 8:[9], 9:[10], 10:[1]})
            sage: s.rauzy_fractal_plot()                                    # not tested (> 1 second)

        #. A non-Pisot reducible substitution with a strange Rauzy fractal::

            sage: s = WordMorphism({1:[3,2], 2:[3,3], 3:[4], 4:[1]})
            sage: s.rauzy_fractal_plot()                                    # not tested (> 1 second)

        #. A substitution with overlapping tiles. We use the options
           ``colormap`` and ``opacity`` to study how the tiles overlap::

            sage: s = WordMorphism('1->213,2->4,3->5,4->1,5->21')
            sage: s.rauzy_fractal_plot()                                    # not tested (> 1 second)
            sage: s.rauzy_fractal_plot(colormap={'1':'red', '4':'purple'})  # not tested (> 1 second)
            sage: s.rauzy_fractal_plot(n=150000,                            # not tested (> 1 second)
            ....:                      opacity={'1':0.1,'2':1,'3':0.1,'4':0.1,'5':0.1})

        #. Funny experiments by playing with the precision of the float numbers used to plot the fractal::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: s.rauzy_fractal_plot(prec=6)                              # not tested
            sage: s.rauzy_fractal_plot(prec=9)                              # not tested
            sage: s.rauzy_fractal_plot(prec=15)                             # not tested
            sage: s.rauzy_fractal_plot(prec=19)                             # not tested
            sage: s.rauzy_fractal_plot(prec=25)                             # not tested

        #. Using the ``translate`` option to plot periodic tilings::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: s.rauzy_fractal_plot(n=10000,                             # not tested (> 1 second)
            ....:                      translate=[(0,0,0),(-1,0,1),(0,-1,1),(1,-1,0),
            ....:                                 (1,0,-1),(0,1,-1),(-1,1,0)])

           ::

            sage: t = WordMorphism("a->aC,b->d,C->de,d->a,e->ab")   # substitution found by Julien Bernat
            sage: V = [vector((0,0,1,0,-1)), vector((0,0,1,-1,0))]                      # needs sage.modules
            sage: S = set(map(tuple, [i*V[0] + j*V[1]                                   # needs sage.modules
            ....:                     for i in [-1,0,1] for j in [-1,0,1]]))
            sage: t.rauzy_fractal_plot(n=10000,                             # not tested (> 1 second)
            ....:                      translate=S, exchange=true)

        #. Using the ``translate`` option to plot arbitrary tilings with the fractal pieces.
           This can be used for example to plot the self-replicating tiling of the Rauzy fractal::

            sage: s = WordMorphism({1:[1,2], 2:[3], 3:[4,3], 4:[5], 5:[6], 6:[1]})
            sage: s.rauzy_fractal_plot()                                    # not tested (> 1 second)
            sage: D = {1: [(0,0,0,0,0,0), (0,1,0,0,0,0)],
            ....:      3: [(0,0,0,0,0,0), (0,1,0,0,0,0)], 6: [(0,1,0,0,0,0)]}
            sage: s.rauzy_fractal_plot(n=30000, translate=D)                # not tested (> 1 second)

        #. Plot the projection of the canonical basis with the fractal::

            sage: s = WordMorphism({1:[2,1], 2:[3], 3:[6,4], 4:[5,1],
            ....:                   5:[6], 6:[7], 7:[8], 8:[9], 9:[1]})
            sage: s.rauzy_fractal_plot(plot_basis=True)                     # not tested (> 1 second)

        TESTS::

            sage: s = WordMorphism('a->ab,b->c,c->d,d->e,e->a')
            sage: s.rauzy_fractal_plot(n=1000, colormap='Set1',                         # needs sage.modules sage.plot
            ....:                      opacity={'a':0.5,'b':1,'c':0.7,'d':0,'e':0.2},
            ....:                      plot_origin=(100,"black"), plot_basis=True,
            ....:                      point_size=2.5)
            Graphics object consisting of 10 graphics primitives

        REFERENCES:

        - [1] Valerie Berthe and Anne Siegel,
          Tilings associated with beta-numeration and substitutions,
          Integers 5 (3), 2005.
          http://www.integers-ejcnt.org/vol5-3.html

        AUTHOR:

            Timo Jolivet (2012-06-16)
        """
        alphabet = self.domain().alphabet()
        size_alphabet = len(alphabet)

        orbit_points = self.rauzy_fractal_points(n=n, exchange=exchange, eig=eig, translate=translate, prec=prec)

        dim_fractal = len(orbit_points[alphabet[0]][0])

        # Manage colors and opacity
        if isinstance(colormap, dict):
            col_dict = colormap

        elif isinstance(colormap, str):
            from matplotlib import cm

            if colormap not in cm.datad:
                raise RuntimeError("color map %s not known (type sorted(colors) for valid names)" % colormap)

            colormap = cm.__dict__[colormap]
            col_dict = {}
            for i, a in enumerate(alphabet):
                col_dict[a] = colormap(float(i) / float(size_alphabet))[:3]

        else:
            raise TypeError("type of option colormap (=%s) must be dict or str" % colormap)

        if opacity is None:
            opacity = {a: 1 for a in alphabet}

        elif not isinstance(opacity, dict):
            raise TypeError("type of option opacity (=%s) must be dict" % opacity)

        # Plot points size
        if point_size is None:
            if dim_fractal == 1 or dim_fractal == 2:
                point_size = 1
            elif dim_fractal == 3:
                point_size = 8

        # Make graphics
        from sage.plot.plot import Graphics
        G = Graphics()

        from sage.plot.point import points

        # 1D plots
        if dim_fractal == 1:
            from sage.plot.plot import plot
            for a in col_dict:
                # We plot only the points with a color in col_dict and with positive opacity
                if (a in col_dict) and (opacity[a] > 0):
                    G += plot([x[0] for x in orbit_points[a]], color=col_dict[a], alpha=opacity[a], thickness=point_size)
            if plot_basis:
                from matplotlib import cm
                from sage.plot.arrow import arrow
                canonical_basis_proj = self.rauzy_fractal_projection(eig=eig, prec=prec)
                for i, a in enumerate(alphabet):
                    x = canonical_basis_proj[a]
                    G += arrow((-1.1, 0), (-1.1, x[0]),
                               color=cm.__dict__["gist_gray"](0.75 * float(i) / float(size_alphabet))[:3])

        # 2D or 3D plots
        else:
            if point_size is None and dim_fractal == 2:
                point_size = 1
            elif point_size is None and dim_fractal == 3:
                point_size = 8

            for a in col_dict:
                # We plot only the points with a color in col_dict and with positive opacity
                if (a in col_dict) and (opacity[a] > 0):
                    G += points(orbit_points[a], color=col_dict[a], alpha=opacity[a], size=point_size)

            if plot_basis:
                from matplotlib import cm
                from sage.plot.arrow import arrow
                canonical_basis_proj = self.rauzy_fractal_projection(eig=eig, prec=prec)
                for i, a in enumerate(alphabet):
                    x = canonical_basis_proj[a]
                    G += arrow([0] * dim_fractal, x,
                               color=cm.__dict__["gist_gray"](0.75 * float(i) / float(size_alphabet))[:3])

        if plot_origin:
            G += points([(0, 0)], size=plot_origin[0], color=plot_origin[1])

        if dim_fractal == 1 or dim_fractal == 2:
            G.set_aspect_ratio(1)

        return G

    def is_growing(self, letter=None):
        r"""
        Return ``True`` if ``letter`` is a growing letter.

        A letter `a` is *growing* for the morphism `s` if the length of the
        iterates of `| s^n(a) |` tend to infinity as `n` goes to infinity.

        INPUT:

        - ``letter`` -- ``None`` or a letter in the domain of ``self``

        .. NOTE::

            If letter is ``None``, this returns ``True`` if ``self`` is
            everywhere growing, i.e., all letters are growing letters (see
            [CassNic10]_), and that ``self`` **must** be an endomorphism.

        EXAMPLES::

            sage: WordMorphism('0->01,1->1').is_growing('0')
            True
            sage: WordMorphism('0->01,1->1').is_growing('1')
            False
            sage: WordMorphism('0->01,1->10').is_growing()
            True
            sage: WordMorphism('0->1,1->2,2->01').is_growing()
            True
            sage: WordMorphism('0->01,1->1').is_growing()
            False

        The domain needs to be equal to the codomain::

            sage: WordMorphism('0->01,1->0,2->1',codomain=Words('012')).is_growing()
            True

        Test of erasing morphisms::

            sage: WordMorphism('0->01,1->').is_growing('0')
            False
            sage: m = WordMorphism('a->bc,b->bcc,c->',codomain=Words('abc'))
            sage: m.is_growing('a')
            False
            sage: m.is_growing('b')
            False
            sage: m.is_growing('c')
            False

        TESTS:

        Make sure that :issue:`31454` is fixed::

            sage: WordMorphism('a->a').is_growing('a')
            False

        REFERENCES:

        ..  [CassNic10] Cassaigne J., Nicolas F. Factor complexity.
            Combinatorics, automata and number theory, 163--247, Encyclopedia
            Math. Appl., 135, Cambridge Univ. Press, Cambridge, 2010.
        """
        if not letter:
            return self.domain().alphabet().cardinality() == len(self.growing_letters())
        else:
            return letter in self.growing_letters()

    def growing_letters(self):
        r"""
        Return the list of growing letters.

        See :meth:`.is_growing` for more information.

        EXAMPLES::

            sage: WordMorphism('0->01,1->10').growing_letters()
            ['0', '1']
            sage: WordMorphism('0->01,1->1').growing_letters()
            ['0']
            sage: WordMorphism('0->01,1->0,2->1',codomain=Words('012')).growing_letters()
            ['0', '1', '2']
            sage: WordMorphism('a->b,b->a').growing_letters()
            []
            sage: WordMorphism('a->b,b->c,c->d,d->c', codomain=Words('abcd')).growing_letters()
            []

        TESTS:

        Make sure that :issue:`31454` is fixed::

            sage: WordMorphism('a->a').growing_letters()
            []
        """
        # Remove letters that vanish, ie sigma^n(letter) is ultimately empty
        immortal = set(self.immortal_letters())
        new_morph = {x: [z for z in self._morph[x] if z in immortal] for x in immortal}

        # Remove cycles of letters
        graph_one = {x: y[0] for x, y in new_morph.items() if len(y) == 1}
        no_loops = set(new_morph)
        for cycle in get_cycles(graph_one.__getitem__, graph_one):
            no_loops.difference_update(cycle)
        new_morph = {x: [z for z in new_morph[x] if z in no_loops] for x in no_loops}

        # NOTE: here we should actually be using the domain made of the
        # remaining letters in new_morph. However, building the corresponding
        # alphabet and finite words cost much more time than using the same
        # domain. Instead we just erase the corresponding letters.
        for a in self._domain.alphabet():
            if a not in new_morph:
                new_morph[a] = self._codomain()

        # Remove letters ending in a cycle
        new_morph = WordMorphism(new_morph, domain=self.domain(), codomain=self.codomain())
        return new_morph.immortal_letters()

    def immortal_letters(self):
        r"""
        Return the list of immortal letters.

        A letter `a` is *immortal* for the morphism `s` if the length of the
        iterates of `| s^n(a) |` is larger than zero as `n` goes to infinity.

        Requires this morphism to be self-composable.

        EXAMPLES::

            sage: WordMorphism('a->a').immortal_letters()
            ['a']
            sage: WordMorphism('a->b,b->a').immortal_letters()
            ['a', 'b']
            sage: WordMorphism('a->abcd,b->cd,c->dd,d->').immortal_letters()
            ['a']
            sage: WordMorphism('a->bc,b->cac,c->de,d->,e->').immortal_letters()
            ['a', 'b']
            sage: WordMorphism('a->', domain=Words('a'), codomain=Words('a')).immortal_letters()
            []

            sage: WordMorphism('a->').immortal_letters()
            []
        """
        if not self.is_self_composable():
            raise TypeError(f'self ({self}) is not a self-composable')

        forward = {}
        backward = {letter: set() for letter in self._morph}
        stack = []
        for letter, image in self._morph.items():
            if not image:
                stack.append(letter)
                forward[letter] = set()
            else:
                simage = set(image)
                forward[letter] = simage
                for occurrence in simage:
                    backward[occurrence].add(letter)

        while stack:
            letter = stack.pop()
            for preimage in backward[letter]:
                forward[preimage].remove(letter)
                if not forward[preimage]:
                    stack.append(preimage)
            del forward[letter]
            del backward[letter]

        return sorted(forward, key=self.domain().alphabet().rank)

    def letter_growth_types(self):
        r"""
        Return the mortal, polynomial and exponential growing letters.

        The growth of `| s^n(a) |` as `n` goes to `\infty` is always of the
        form `\alpha^n n^\beta` (where `\alpha` is a Perron number and
        `\beta` an integer).

        Without doing any linear algebra three cases can be differentiated:
        mortal (ultimately empty or `\alpha=0`); polynomial (`\alpha=1`);
        exponential (`\alpha > 1`). This is what is done in this method.

        It requires this morphism to be an endomorphism.

        OUTPUT:

        The output is a 3-tuple of lists (mortal, polynomial, exponential)
        where:

        - ``mortal`` -- list of mortal letters
        - ``polynomial`` -- list of lists where ``polynomial[i]`` is the
          list of letters with growth `n^i`
        - ``exponential`` -- list of at least exponentionally growing letters

        EXAMPLES::

            sage: s = WordMorphism('a->abc,b->bc,c->c')
            sage: mortal, poly, expo = s.letter_growth_types()
            sage: mortal
            []
            sage: poly
            [['c'], ['b'], ['a']]
            sage: expo
            []

        When three mortal letters (c, d, and e), and two letters (a, b) are
        not growing::

            sage: s = WordMorphism('a->bc,b->cac,c->de,d->,e->')
            sage: s^20
            WordMorphism: a->cacde, b->debcde, c->, d->, e->
            sage: mortal, poly, expo = s.letter_growth_types()
            sage: mortal
            ['c', 'd', 'e']
            sage: poly
            [['a', 'b']]
            sage: expo
            []

        ::

            sage: s = WordMorphism('a->abcd,b->bc,c->c,d->a')
            sage: mortal, poly, expo = s.letter_growth_types()
            sage: mortal
            []
            sage: poly
            [['c'], ['b']]
            sage: expo
            ['a', 'd']

        TESTS::

            sage: s = WordMorphism('a->a')
            sage: s.letter_growth_types()
            ([], [['a']], [])

        ::

            sage: s = WordMorphism('a->b,b->a')
            sage: s.letter_growth_types()
            ([], [['a', 'b']], [])

        ::

            sage: s = WordMorphism('a->abcd,b->cd,c->dd,d->')
            sage: s.letter_growth_types()
            (['b', 'c', 'd'], [['a']], [])

        ::

            sage: s = WordMorphism('a->', domain=Words('a'), codomain=Words('a'))
            sage: s.letter_growth_types()
            (['a'], [], [])
        """
        immortal = set(self.immortal_letters())
        mortal = [a for a in self.domain().alphabet()
                  if a not in immortal]

        # Starting with degree d=0, search for letters with polynomial
        # growth of degree d.
        polynomial = []
        m = {a: [b for b in self.image(a) if b in immortal] for a in immortal}
        while True:
            # Construct the permutation of letters containing all letters whose
            # iterated images under morphism m is always of length 1.
            not_growing = {a: image_a[0] for a, image_a in m.items()
                           if len(image_a) == 1}
            preimages = {}
            roots = []
            for k, v in not_growing.items():
                if v not in not_growing:
                    roots.append(v)
                if v not in preimages:
                    preimages[v] = []
                preimages[v].append(k)

            while roots:
                v = roots.pop()
                for k in preimages.get(v):
                    del not_growing[k]
                    if k in preimages:
                        roots.append(k)

            # The letters inside not_growing are the ones with polynomial
            # growth d. If there is none, then the remaining letters in m
            # have exponential growth.
            if not not_growing:
                break
            polynomial.append(list(not_growing))

            # clean the morphism m for the next iteration by removing the
            # letters with polynomial growth degree d
            m = {a: [b for b in L if b not in not_growing] for a, L in m.items()
                 if a not in not_growing}

        exponential = list(m)

        # sort the letters as in the input alphabet if possible
        A = self.domain().alphabet()
        try:
            rank = A.rank
        except AttributeError:
            pass
        else:
            mortal.sort(key=rank)
            for letters in polynomial:
                letters.sort(key=rank)
            exponential.sort(key=rank)

        return mortal, polynomial, exponential

    def abelian_rotation_subspace(self):
        r"""
        Return the subspace on which the incidence matrix of ``self`` acts by
        roots of unity.

        EXAMPLES::

            sage: # needs sage.modules
            sage: WordMorphism('0->1,1->0').abelian_rotation_subspace()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            sage: WordMorphism('0->01,1->10').abelian_rotation_subspace()
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: WordMorphism('0->01,1->1').abelian_rotation_subspace()
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
            sage: WordMorphism('1->122,2->211').abelian_rotation_subspace()
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -1]
            sage: WordMorphism('0->1,1->102,2->3,3->4,4->2').abelian_rotation_subspace()
            Vector space of degree 5 and dimension 3 over Rational Field
            Basis matrix:
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]

        The domain needs to be equal to the codomain::

            sage: WordMorphism('0->1,1->',codomain=Words('01')).abelian_rotation_subspace()         # needs sage.modules
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        if not self.domain() == self.codomain():
            raise TypeError("self (=%s) is not an endomorphism" % self)

        if self.domain().alphabet().cardinality() == Infinity:
            raise ValueError("the alphabet is infinite")

        M = self.incidence_matrix()
        p = M.charpoly().factor()
        basis = []
        for factor in p:
            if factor[0].is_cyclotomic():
                basis.extend((factor[0])(M).right_kernel().basis())

        return M.column_ambient_module(base_ring=QQ).subspace(basis)

    def is_injective(self):
        """
        Return whether this morphism is injective.

        ALGORITHM:

        Uses a version of :wikipedia:`Sardinas–Patterson_algorithm`.
        Time complexity is on average quadratic with regards to the size of the
        morphism.

        EXAMPLES::

            sage: WordMorphism('a->0,b->10,c->110,d->111').is_injective()
            True
            sage: WordMorphism('a->00,b->01,c->012,d->20001').is_injective()
            False
        """
        def check(u, v):
            if u.is_prefix(v):
                tail = v[u.length():]
                if tail not in tails:
                    tails.add(tail)
                    todo.append(tail)

        if self.is_erasing():
            return False
        images = self.images()
        tails = set()
        todo = []

        for i in range(len(images)):
            for j in range(i + 1, len(images)):
                if images[i] == images[j]:
                    return False
                check(images[i], images[j])
                check(images[j], images[i])
        while todo:
            u = todo.pop()
            for v in images:
                if u == v:
                    return False
                check(u, v)
                check(v, u)

        return True

    def is_pushy(self, w=None):
        r"""
        Return whether the language `\{m^n(w) | n \ge 0\}` is pushy,
        where `m` is this morphism and `w` is a word inputted as a parameter.

        Requires this morphism to be an endomorphism.

        A language created by iterating a morphism is pushy, if its words
        contain an infinite number of factors containing no growing letters. It
        turns out that this is equivalent to having at least one infinite
        repetition containing no growing letters.

        See :meth:`infinite_repetitions_primitive_roots` and :meth:`is_growing`.

        INPUT:

        - ``w`` -- finite iterable (default: ``self.domain().alphabet()``);
          represents a word used to start the language

        EXAMPLES::

            sage: WordMorphism('a->abca,b->bc,c->').is_pushy()
            False
            sage: WordMorphism('a->abc,b->,c->bcb').is_pushy()
            True
        """
        return bool(self.infinite_repetitions_primitive_roots(w, False))

    def is_unboundedly_repetitive(self, w=None):
        r"""
        Return whether the language `\{m^n(w) | n \ge 0\}` is unboundedly repetitive,
        where `m` is this morphism and `w` is a word inputted as a parameter.

        Requires this morphism to be an endomorphism.

        A language created by iterating a morphism is unboundedly repetitive, if
        it has at least one infinite repetition containing at least one growing
        letter.

        See :meth:`infinite_repetitions_primitive_roots` and :meth:`is_growing`.

        INPUT:

        - ``w`` -- finite iterable (default: ``self.domain().alphabet()``);
          represents a word used to start the language

        EXAMPLES::

            sage: WordMorphism('a->abca,b->bc,c->').is_unboundedly_repetitive()
            True
            sage: WordMorphism('a->abc,b->,c->bcb').is_unboundedly_repetitive()
            False
        """
        return bool(self.infinite_repetitions_primitive_roots(w, True))

    def is_repetitive(self, w=None):
        r"""
        Return whether the language `\{m^n(w) | n \ge 0\}` is repetitive,
        where `m` is this morphism and `w` is a word inputted as a parameter.

        Requires this morphism to be an endomorphism.

        A language is repetitive, if for each positive integer `k` there exists
        a word `u` such that `u^k` is a factor of some word of the language.

        It turns out that for languages created by iterating a morphism this is
        equivalent to having at least one infinite repetition (this property is
        also known as strong repetitiveness).

        See :meth:`infinite_repetitions_primitive_roots`.

        INPUT:

        - ``w`` -- finite iterable (default: ``self.domain().alphabet()``);
          represents a word used to start the language

        EXAMPLES:

        This method can be used to check whether a purely morphic word is not
        k-power free for all positive integers k. For example, the language
        containing just the Thue-Morse word and its prefixes is not repetitive,
        since the Thue-Morse word is cube-free::

            sage: WordMorphism('a->ab,b->ba').is_repetitive('a')
            False

        Similarly, the Hanoi word is square-free::

            sage: WordMorphism('a->aC,A->ac,b->cB,B->cb,c->bA,C->ba').is_repetitive('a')
            False

        However, this method solves a more general problem, as it can be called
        on any morphism `m` and with any word `w`::

            sage: WordMorphism('a->c,b->cda,c->a,d->abc').is_repetitive('bd')
            True
        """
        return self.is_pushy(w) or self.is_unboundedly_repetitive(w)

    def infinite_repetitions_primitive_roots(self, w=None, allow_growing=None):
        r"""
        Return the set of primitive roots (up to conjugacy) of infinite
        repetitions from the language `\{m^n(w) | n \ge 0\}`, where `m` is this
        morphism and `w` is a word inputted as a parameter.

        Requires this morphism to be an endomorphism.

        The word `v^\omega` is an infinite repetition (in other words, an
        infinite periodic factor) of a language, if `v` is a non-empty word and
        for each positive integer `k` the word `v^k` is a factor of some word
        from the language. It turns out that a language created by iterating a
        morphism has a finite number of primitive roots of infinite repetitions.

        If `v` is a primitive root of an infinite repetition, then all its
        conjugations are also primitive roots of an infinite repetition. For
        simplicity's sake this method returns only the lexicographically minimal
        one from each conjugacy class.

        INPUT:

        - ``w`` -- finite iterable (default: ``self.domain().alphabet()``);
          represents a word used to start the language

        - ``allow_growing`` -- boolean or ``None`` (default: ``None``); if
          ``False``, return only the primitive roots that contain no growing
          letters. If ``True``, return only the primitive roots that contain at
          least one growing letter. If ``None``, return both.

        ALGORITHM:

        The algorithm used is described in detail in [KS2015]_.

        EXAMPLES::

            sage: m = WordMorphism('a->aba,b->aba,c->cd,d->e,e->d')
            sage: inf_reps = m.infinite_repetitions_primitive_roots('ac')
            sage: sorted(inf_reps)
            [word: aab, word: de]

        ``allow_growing`` parameter::

            sage: sorted(m.infinite_repetitions_primitive_roots('ac', True))
            [word: aab]
            sage: sorted(m.infinite_repetitions_primitive_roots('ac', False))
            [word: de]

        Incomplete check that these words are indeed the primitive roots of
        infinite repetitions::

            sage: SL = m._language_naive(10, Word('ac'))
            sage: all(x in SL for x in inf_reps)
            True
            sage: all(x^2 in SL for x in inf_reps)
            True
            sage: all(x^3 in SL for x in inf_reps)
            True

        Large example::

            sage: m = WordMorphism('a->1b5,b->fcg,c->dae,d->432,e->678,f->f,g->g,1->2,2->3,3->4,4->1,5->6,6->7,7->8,8->5')
            sage: sorted(m.infinite_repetitions_primitive_roots('a'))
            [word: 1432f2143f3214f4321f, word: 5678g8567g7856g6785g]

        TESTS::

            sage: m = WordMorphism('a->Cab,b->1c1,c->E2bd5,d->BbaA,5->6,6->7,7->8,8->9,9->5,1->2,2->1,A->B,B->C,C->D,D->E,E->')
            sage: sorted(m.infinite_repetitions_primitive_roots())
            [word: 1, word: 1519181716, word: 2, word: 2529282726]

            sage: m = WordMorphism('a->b,b->b', codomain=FiniteWords('ab'))
            sage: m.infinite_repetitions_primitive_roots()
            set()

            sage: m = WordMorphism('c->d,d->c,e->fc,f->ed')
            sage: sorted(m.infinite_repetitions_primitive_roots())
            [word: c, word: d]

            sage: m = WordMorphism('a->bcb,b->ada,c->d,d->c')
            sage: sorted(m.infinite_repetitions_primitive_roots())
            [word: ad, word: bc]

            sage: m = WordMorphism('b->c,c->bcb')
            sage: sorted(m.infinite_repetitions_primitive_roots())
            [word: bc]

            sage: m = WordMorphism('a->abc,b->dab,c->abc,d->dab')
            sage: sorted(m.infinite_repetitions_primitive_roots())
            [word: ababcd]
        """
        def impl_no_growing(g, k):
            U = {}
            for x in unbounded:
                xg = g.image(x)
                for i, y in enumerate(reversed(xg)):
                    if y in unbounded:
                        break
                U[x] = y, xg[xg.length() - i:]
            for cycle in get_cycles(lambda x: U[x][0], domain=unbounded):
                if all(not U[x][1] for x in cycle):
                    continue
                gq = gb**len(cycle)
                for cyc in g.domain()(cycle).conjugates_iterator():
                    u = g.domain()()
                    for x in cyc:
                        u = U[x][1] + gb(u)
                    inf_rep = g.domain()()
                    history = set()
                    while u not in history:
                        history.add(u)
                        inf_rep += u
                        u = gq(u)
                    yield k(inf_rep.primitive()).primitive()

        if w is None:
            w = self._morph
        reach = self._language_naive(2, self._domain(w))
        f = self.restrict_domain([x[0] for x in reach])
        f._codomain = f._domain
        g, _, k, _ = f.simplify_until_injective()
        g._codomain = g._domain
        unbounded = set(g.growing_letters())
        result = set()

        if allow_growing is not True:
            gb = g.restrict_domain(set(g._morph) - unbounded)
            for x in impl_no_growing(g, k):  # UR.
                result.add(x.minimal_conjugate())
            for x in impl_no_growing(g.reversal(), k.reversal()):  # UL.
                result.add(self.domain()(reversed(x)).minimal_conjugate())

        if allow_growing is not False:
            for periodic_orbit in g.periodic_points():
                gq = g**len(periodic_orbit)
                for periodic_point in periodic_orbit:
                    # Check if this periodic point is a periodic infinite word.
                    periodic_point = periodic_point[:1]
                    occurred = set(periodic_point)
                    one_unbounded_twice = False
                    for _ in g.domain().alphabet():
                        previous_length = periodic_point.length()
                        periodic_point = gq(periodic_point)
                        for i, letter in enumerate(periodic_point[previous_length:]):
                            if letter in unbounded:
                                if letter in occurred:
                                    one_unbounded_twice = True
                                    break
                                occurred.add(letter)
                        if one_unbounded_twice:
                            break
                    if not one_unbounded_twice or letter != periodic_point[0]:
                        break
                    v = periodic_point[:previous_length + i]
                    vq = gq(v)
                    m = 0
                    while vq[m * v.length(): (m + 1) * v.length()] == v:
                        m += 1
                    if m * v.length() != vq.length():
                        break
                    result.add(k(v).primitive().minimal_conjugate())

        return result

    def simplify_alphabet_size(self, Z=None):
        r"""
        If this morphism is simplifiable, return morphisms `h` and `k` such that
        this morphism is simplifiable with respect to `h` and `k`, otherwise
        raise :exc:`ValueError`.

        This method is quite fast if this morphism is non-injective, but very
        slow if it is injective.

        Let `f: X^* \rightarrow Y^*` be a morphism. Then `f` is simplifiable
        with respect to morphisms `h: X^* \rightarrow Z^*` and
        `k: Z^* \rightarrow Y^*`, if `f = k \circ h` and `|Z| < |X|`. If also
        `Y \subseteq X`, then the morphism `g: Z^* \rightarrow Z^* = h \circ k`
        is a simplification of `f` (with respect to `h` and `k`).

        Loosely speaking, a morphism is simplifiable if it contains "more letters
        than is needed". Non-injectivity implies simplifiability. Simplification
        preserves some properties of the original morphism (e.g. repetitiveness).

        For more information see Section 3 in [KO2000]_.

        INPUT:

        - ``Z`` -- iterable (default: ``self.domain().alphabet()``) whose
          elements are used as an alphabet for the simplification

        EXAMPLES:

        Example of a simplifiable (non-injective) morphism::

            sage: f = WordMorphism('a->aca,b->badc,c->acab,d->adc')
            sage: h, k = f.simplify_alphabet_size('xyz'); h, k
            (WordMorphism: a->x, b->zy, c->xz, d->y, WordMorphism: x->aca, y->adc, z->b)
            sage: k * h == f
            True
            sage: g = h * k; g
            WordMorphism: x->xxzx, y->xyxz, z->zy

        Example of a simplifiable (injective) morphism::

            sage: f = WordMorphism('a->abcc,b->abcd,c->abdc,d->abdd')
            sage: h, k = f.simplify_alphabet_size('xyz'); h, k
            (WordMorphism: a->xyy, b->xyz, c->xzy, d->xzz, WordMorphism: x->ab, y->c, z->d)
            sage: k * h == f
            True
            sage: g = h * k; g
            WordMorphism: x->xyyxyz, y->xzy, z->xzz

        Example of a non-simplifiable morphism::

            sage: WordMorphism('a->aa').simplify_alphabet_size()
            Traceback (most recent call last):
            ...
            ValueError: self (a->aa) is not simplifiable

        Example of an erasing morphism::

            sage: f = WordMorphism('a->abc,b->cc,c->')
            sage: h, k = f.simplify_alphabet_size(); h, k
            (WordMorphism: a->a, b->b, c->, WordMorphism: a->abc, b->cc)
            sage: k * h == f
            True
            sage: g = h * k; g
            WordMorphism: a->ab, b->

        Example of a morphism, that is not an endomorphism::

            sage: f = WordMorphism('a->xx,b->xy,c->yx,d->yy')
            sage: h, k = f.simplify_alphabet_size(NN); h, k
            (WordMorphism: a->00, b->01, c->10, d->11, WordMorphism: 0->x, 1->y)
            sage: k * h == f
            True
            sage: len(k.domain().alphabet()) < len(f.domain().alphabet())
            True
        """
        def try_create_h(f, k):
            h = {}
            for letter1, image1 in f.items():
                image3 = []
                while image1:
                    for letter2, image2 in k.items():
                        if image2.is_prefix(image1):
                            image1 = image1[image2.length():]
                            image3.append(letter2)
                            break
                    else:  # nobreak
                        return None
                h[letter1] = image3
            return h

        X = self.domain().alphabet()
        Y = self.codomain().alphabet()
        f = self._morph

        if self.is_erasing():  # Trivial case #1.
            k = {letter: image for letter, image in f.items() if image}
            h = {letter: [letter] if image else [] for letter, image in f.items()}
        elif len(Y) < len(X):  # Trivial case #2.
            k = {x: [y] for x, y in zip(X, Y)}
            k_inverse = dict(zip(Y, X))
            h = {x: [k_inverse[y] for y in image] for x, image in f.items()}
        elif not self.is_injective():  # Non-trivial but a fast case.
            k = dict(f)
            to_do = set(k)
            while to_do:
                to_remove = []
                # min() and remove() instead of pop() to have deterministic output.
                letter1 = min(to_do)
                to_do.remove(letter1)
                image1 = k[letter1]
                for letter2, image2 in k.items():
                    if letter1 == letter2:
                        continue
                    if image1 == image2:
                        to_remove.append(letter2)
                        to_do.discard(letter2)
                    elif image1.is_prefix(image2):
                        k[letter2] = image2[image1.length():]
                        to_do.add(letter2)
                    elif image2.is_prefix(image1):
                        k[letter1] = image1[image2.length():]
                        to_do.add(letter1)
                        break
                for letter in to_remove:
                    del k[letter]
            h = try_create_h(f, k)
        else:  # Non-trivial and a slow case.
            factors = set()
            for image in f.values():
                factors.update(x.primitive() for x in image.factor_iterator())
            factors.remove(self.codomain()())
            factors = sorted(factors)  # For deterministic output.
            from itertools import combinations
            for comb in combinations(factors, len(X) - 1):
                if any(x.is_proper_prefix(y) for x in comb for y in comb):
                    continue
                k = dict(zip(X, comb))
                h = try_create_h(f, k)
                if h:
                    break
            else:  # nobreak
                raise ValueError(f'self ({self}) is not simplifiable')

        k = WordMorphism(k, codomain=self.codomain())
        h = WordMorphism(h, domain=self.domain(), codomain=k.domain())

        if Z is not None:  # Custom alphabet.
            old_Z_star = k.domain()
            old_Z = old_Z_star.alphabet()
            Z = [z for z, _ in zip(Z, old_Z)]
            if len(Z) < len(old_Z):
                raise ValueError(f'Z should have length at least {len(old_Z)}, is {len(Z)}')
            Z_star = FiniteWords(Z)
            h_new = {old: [new] for old, new in zip(old_Z, Z)}
            k_new = {new: [old] for new, old in zip(Z, old_Z)}
            h_new = WordMorphism(h_new, domain=old_Z_star, codomain=Z_star)
            k_new = WordMorphism(k_new, domain=Z_star, codomain=old_Z_star)
            h = h_new * h
            k = k * k_new

        return h, k

    def simplify_until_injective(self):
        r"""
        Return a quadruplet `(g, h, k, i)`, where `g` is an injective
        simplification of this morphism with respect to `h`, `k` and `i`.

        Requires this morphism to be an endomorphism.

        This methods basically calls :meth:`simplify_alphabet_size` until the
        returned simplification is injective. If this morphism is already
        injective, a quadruplet `(g, h, k, i)` is still returned, where `g`
        is this morphism, `h` and `k` are the identity morphisms and `i` is 0.

        Let `f: X^* \rightarrow Y^*` be a morphism and `Y \subseteq X`. Then
        `g: Z^* \rightarrow Z^*` is an injective simplification of `f` with
        respect to morphisms `h: X^* \rightarrow Z^*` and
        `k: Z^* \rightarrow Y^*` and a positive integer `i`, if `g` is
        injective, `|Z| < |X|`, `g^i = h \circ k` and `f^i = k \circ h`.

        For more information see Section 4 in [KO2000]_.

        EXAMPLES::

            sage: f = WordMorphism('a->abc,b->a,c->bc')
            sage: g, h, k, i = f.simplify_until_injective(); g, h, k, i
            (WordMorphism: a->aa, WordMorphism: a->aa, b->a, c->a, WordMorphism: a->abc, 2)
            sage: g.is_injective()
            True
            sage: g**i == h * k
            True
            sage: f**i == k * h
            True
        """
        if not self.is_endomorphism():
            raise TypeError(f'self ({self}) is not an endomorphism')

        g = self
        h = self.domain().identity_morphism()
        k = self.codomain().identity_morphism()
        i = 0
        while not g.is_injective():
            h_new, k_new = g.simplify_alphabet_size()
            g, h, k, i = h_new * k_new, h_new * h, k * k_new, i + 1
        return g, h, k, i
