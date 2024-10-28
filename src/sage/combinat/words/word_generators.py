r"""
Common words

AUTHORS:

- Franco Saliola (2008-12-17): merged into sage
- Sébastien Labbé (2008-12-17): merged into sage
- Arnaud Bergeron (2008-12-17): merged into sage
- Amy Glen (2008-12-17): merged into sage
- Sébastien Labbé (2009-12-19): Added S-adic words (:issue:`7543`)

USE:

To see a list of all word constructors, type ``words.`` and then press
the :kbd:`Tab` key. The documentation for each constructor includes
information about each word, which provides a useful reference.

REFERENCES:

.. [AC03] \B. Adamczewski, J. Cassaigne, On the transcendence of real
   numbers with a regular expansion, J. Number Theory 103 (2003)
   27--37.

.. [BmBGL07] \A. Blondin-Massé, S. Brlek, A. Glen, and S. Labbé. On the
   critical exponent of generalized Thue-Morse words. *Discrete Math.
   Theor. Comput.  Sci.* 9 (1):293--304, 2007.

.. [BmBGL09] \A. Blondin-Massé, S. Brlek, A. Garon, and S. Labbé. Christoffel
   and Fibonacci Tiles, DGCI 2009, Montreal, to appear in LNCS.

.. [Loth02] \M. Lothaire, Algebraic Combinatorics On Words, vol. 90 of
   Encyclopedia of Mathematics and its Applications, Cambridge
   University Press, U.K., 2002.

.. [Fogg] Pytheas Fogg,
   https://www.lirmm.fr/arith/wiki/PytheasFogg/S-adiques.

EXAMPLES::

    sage: t = words.ThueMorseWord(); t
    word: 0110100110010110100101100110100110010110...
"""
# ****************************************************************************
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>,
#                          Sebastien Labbe <slabqc@gmail.com>,
#                          Arnaud Bergeron <abergeron@gmail.com>,
#                          Amy Glen <amy.glen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import Iterable
from itertools import cycle, count
from random import randint
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.combinat.words.abstract_word import Word_class
from sage.combinat.words.word import FiniteWord_list
from sage.combinat.words.finite_word import FiniteWord_class, Factorization
from sage.combinat.words.words import FiniteWords, InfiniteWords
from sage.combinat.words.morphism import WordMorphism
from sage.arith.misc import gcd
from sage.misc.decorators import rename_keyword
from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.real_mpfr', 'RR')


def _build_tab(sym, tab, W):
    r"""
    Internal function building a coding table for the ``phi_inv_tab`` function.

    TESTS::

        sage: from sage.combinat.words.word_generators import _build_tab
        sage: _build_tab(1, [], Words([1, 2]))
        [1]
        sage: _build_tab(1, [1], Words([1, 2]))
        [1, 2]
        sage: _build_tab(2, [1], Words([1, 2]))
        [2, 2]
        sage: _build_tab(2, [1, 2], Words([1, 2]))
        [2, 2, 1]
        sage: _build_tab(1, [2, 2], Words([1, 2]))
        [1, 1, 2]
    """
    c = W.alphabet().cardinality()
    res = [sym]
    if len(tab) == 0:
        return res
    if sym == 1:
        res += tab
        res[1] = (res[1] % c) + 1
        return res
    w = W([sym]).delta_inv(W, tab[0])
    w = w[1:]
    res.append((w[-1] % c) + 1)
    for i in range(1, len(tab)):
        w = w.delta_inv(W, tab[i])
        res.append((w[-1] % c) + 1)
    return res


class LowerChristoffelWord(FiniteWord_list):
    r"""
    Return the lower Christoffel word of slope `p/q`, where `p` and
    `q` are relatively prime nonnegative integers, over the given
    two-letter alphabet.

    The *Christoffel word of slope `p/q`* is obtained from the
    Cayley graph of `\ZZ/(p+q)\ZZ` with generator `q` as
    follows. If `u \rightarrow v` is an edge in the Cayley graph, then
    `v = u + p \mod{p+q}`. Label the edge `u \rightarrow v` by
    ``alphabet[1]`` if `u < v` and ``alphabet[0]`` otherwise. The Christoffel
    word is the word obtained by reading the edge labels along the
    cycle beginning from 0.

    EXAMPLES::

        sage: words.LowerChristoffelWord(4,7)
        word: 00100100101

    ::

        sage: words.LowerChristoffelWord(4,7,alphabet='ab')
        word: aabaabaabab

    TESTS::

        sage: words.LowerChristoffelWord(1,0)
        word: 1
        sage: words.LowerChristoffelWord(0,1,'xy')
        word: x
        sage: words.LowerChristoffelWord(1,1)
        word: 01
    """

    def __init__(self, p, q, alphabet=(0, 1), algorithm='cf'):
        r"""
        INPUT:

        - ``p`` -- integer coprime with `q`
        - ``q`` -- integer coprime with `p`
        - ``alphabet`` -- sequence of two elements (default: (0, 1))
        - ``algorithm`` -- construction method (default: ``'cf'``).
          It can be one of the following:

          - ``'linear'`` -- linear algorithm in the length of the word
          - ``'cf'`` -- fast method using continued fraction

        TESTS::

            sage: words.ChristoffelWord(9, 4, algorithm='linear')
            word: 0110110110111
            sage: words.ChristoffelWord(9, 4, algorithm='cf')
            word: 0110110110111
            sage: words.ChristoffelWord(4, 9, algorithm='linear')
            word: 0001001001001
            sage: words.ChristoffelWord(4, 9, algorithm='cf')
            word: 0001001001001

        ::

            sage: words.LowerChristoffelWord(4,8)
            Traceback (most recent call last):
            ...
            ValueError: 4 and 8 are not relatively prime
            sage: words.LowerChristoffelWord(17, 39, 'xyz')
            Traceback (most recent call last):
            ...
            ValueError: alphabet must contain exactly two distinct elements
            sage: w = words.LowerChristoffelWord(4,7)
            sage: w2 = loads(dumps(w))
            sage: w == w2
            True
            sage: type(w2)
            <class 'sage.combinat.words.word_generators.LowerChristoffelWord'>
            sage: _ = w2.standard_factorization() # hackish test for self.__p and self.__q
        """
        if len(set(alphabet)) != 2:
            raise ValueError("alphabet must contain exactly two distinct elements")
        # Compute gcd of p, q; raise TypeError if not 1.
        if gcd(p, q) != 1:
            raise ValueError("%s and %s are not relatively prime" % (p, q))
        # Compute the Christoffel word
        if algorithm == 'linear':
            w = []
            u = 0
            if (p, q) == (0, 1):
                w = [alphabet[0]]
            else:
                for i in range(p + q):
                    v = (u + p) % (p + q)
                    new_letter = alphabet[0] if u < v else alphabet[1]
                    w.append(new_letter)
                    u = v
        elif algorithm == 'cf':
            if (p, q) == (0, 1):
                w = [alphabet[0]]
            elif (p, q) == (1, 0):
                w = [alphabet[1]]
            else:
                from sage.rings.rational_field import QQ
                cf = QQ((p, q)).continued_fraction_list()
                u = [alphabet[0]]
                v = [alphabet[1]]
                # do not consider the first zero if p < q
                start = 1 if p < q else 0
                for i in range(start, len(cf)-1):
                    if i % 2 == 0:
                        u = u + v * cf[i]
                    else:
                        v = u * cf[i] + v
                i = len(cf)-1
                if i % 2 == 0:
                    u = u + v * (cf[i]-1)
                else:
                    v = u * (cf[i]-1) + v
                w = u + v
        else:
            raise ValueError('Unknown algorithm (=%s)' % algorithm)
        super().__init__(FiniteWords(alphabet), w)
        self.__p = p
        self.__q = q

    def markoff_number(self):
        r"""
        Return the Markoff number associated to the Christoffel word ``self``.

        The *Markoff number* of a Christoffel word `w` is `trace(M(w))/3`,
        where `M(w)` is the `2\times 2` matrix obtained by applying the
        morphism:
        0 -> matrix(2,[2,1,1,1])
        1 -> matrix(2,[5,2,2,1])

        EXAMPLES::

            sage: w0 = words.LowerChristoffelWord(4,7)
            sage: w1, w2 = w0.standard_factorization()
            sage: (m0,m1,m2) = (w.markoff_number() for w in (w0,w1,w2))                 # needs sage.modules
            sage: (m0,m1,m2)                                                            # needs sage.modules
            (294685, 13, 7561)
            sage: m0**2 + m1**2 + m2**2 == 3*m0*m1*m2                                   # needs sage.modules
            True
        """
        from sage.matrix.constructor import matrix
        eta = {0:matrix(2,[2,1,1,1]), 1:matrix(2,[5,2,2,1])}
        M = matrix(2,[1,0,0,1])
        for a in self:
            M *= eta[a]
        return M.trace()/3

    def standard_factorization(self):
        r"""
        Return the standard factorization of the Christoffel word ``self``.

        The *standard factorization* of a Christoffel word `w` is the
        unique factorization of `w` into two Christoffel words.

        EXAMPLES::

            sage: w = words.LowerChristoffelWord(5,9)
            sage: w
            word: 00100100100101
            sage: w1, w2 = w.standard_factorization()
            sage: w1
            word: 001
            sage: w2
            word: 00100100101

        ::

            sage: w = words.LowerChristoffelWord(51,37)
            sage: w1, w2 = w.standard_factorization()
            sage: w1
            word: 0101011010101101011
            sage: w2
            word: 0101011010101101011010101101010110101101...
            sage: w1 * w2 == w
            True
        """
        p, q = self.__p, self.__q
        index = 0
        u = 0
        for i in range(p + q):
            v = (u+p) % (p+q)
            if v == 1:
                index = i
                break
            u = v
        w1, w2 = self[:index+1], self[index+1:]
        w10 = w1.number_of_letter_occurrences(0)
        w11 = w1.number_of_letter_occurrences(1)
        w20 = w2.number_of_letter_occurrences(0)
        w21 = w2.number_of_letter_occurrences(1)
        return Factorization([LowerChristoffelWord(w11,w10),
                              LowerChristoffelWord(w21,w20)])

    def __reduce__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.word_generators import LowerChristoffelWord
            sage: w = LowerChristoffelWord(5,7)
            sage: w.__reduce__()
            (<class 'sage.combinat.words.word_generators.LowerChristoffelWord'>, (5, 7, {0, 1}))
        """
        return self.__class__, (self.__p, self.__q, self.parent().alphabet())


class WordGenerator:
    r"""
    Constructor of several famous words.

    EXAMPLES::

        sage: words.ThueMorseWord()
        word: 0110100110010110100101100110100110010110...

    ::

        sage: words.FibonacciWord()
        word: 0100101001001010010100100101001001010010...

    ::

        sage: words.ChristoffelWord(5, 8)
        word: 0010010100101

    ::

        sage: words.RandomWord(10, 4)    # not tested random
        word: 1311131221

    ::

        sage: words.CodingOfRotationWord(alpha=0.618, beta=0.618)
        word: 1010110101101101011010110110101101101011...

    ::

        sage: tm = WordMorphism('a->ab,b->ba')
        sage: fib = WordMorphism('a->ab,b->a')
        sage: tmword = words.ThueMorseWord([0, 1])
        sage: from itertools import repeat
        sage: words.s_adic(tmword, repeat('a'), {0:tm, 1:fib})
        word: abbaababbaabbaabbaababbaababbaabbaababba...

    .. NOTE::

        To see a list of all word constructors, type ``words.`` and then
        hit the :kbd:`Tab` key. The documentation for each constructor
        includes information about each word, which provides a useful
        reference.

    TESTS::

        sage: from sage.combinat.words.word_generators import WordGenerator
        sage: words2 = WordGenerator()
        sage: type(loads(dumps(words2)))
        <class 'sage.combinat.words.word_generators.WordGenerator'>
    """
    def ThueMorseWord(self, alphabet=(0, 1), base=2):
        r"""
        Return the (Generalized) Thue-Morse word over the given alphabet.

        There are several ways to define the Thue-Morse word `t`.
        We use the following definition: `t[n]` is the sum modulo `m` of
        the digits in the given base expansion of `n`.

        See [BmBGL07]_, [Brlek89]_, and [MH38]_.

        INPUT:

        - ``alphabet`` -- (default: ``(0, 1)``) any container that is suitable
          to build an instance of OrderedAlphabet (``list``, ``tuple``, ``str``, ...)

        - ``base`` -- integer (default: 2); greater than or equal to 2

        EXAMPLES:

        Thue-Morse word::

            sage: t = words.ThueMorseWord(); t
            word: 0110100110010110100101100110100110010110...

        Thue-Morse word on other alphabets::

            sage: t = words.ThueMorseWord('ab'); t
            word: abbabaabbaababbabaababbaabbabaabbaababba...

        ::

            sage: t = words.ThueMorseWord(['L1', 'L2'])
            sage: t[:8]
            word: L1,L2,L2,L1,L2,L1,L1,L2

        Generalized Thue Morse word::

            sage: words.ThueMorseWord(alphabet=(0,1,2), base=2)
            word: 0112122012202001122020012001011212202001...
            sage: t = words.ThueMorseWord(alphabet=(0,1,2), base=5); t
            word: 0120112012201200120112012120122012001201...
            sage: t[100:130].critical_exponent()
            10/3

        TESTS::

            sage: words.ThueMorseWord(alphabet='ab', base=1)
            Traceback (most recent call last):
            ...
            ValueError: base (=1) and len(alphabet) (=2) must be at least 2

        REFERENCES:

        .. [Brlek89] Brlek, S. 1989. «Enumeration of the factors in the Thue-Morse
           word», *Discrete Appl. Math.*, vol. 24, p. 83--96.

        .. [MH38] Morse, M., et G. A. Hedlund. 1938. «Symbolic dynamics»,
           *American Journal of Mathematics*, vol. 60, p. 815--866.
        """
        W = InfiniteWords(alphabet)
        alphabet = W.alphabet()
        m = alphabet.cardinality()
        if base < 2 or m < 2:
            raise ValueError("base (=%s) and len(alphabet) (=%s) must be at least 2" % (base, m))
        from functools import partial
        f = partial(self._ThueMorseWord_nth_digit, alphabet=alphabet, base=base)
        return W(f, datatype='callable')

    def _ThueMorseWord_nth_digit(self, n, alphabet=(0, 1), base=2):
        r"""
        Return the `n`-th letter of the (Generalized) Thue-Morse word.

        The `n`-th digit of the Thue-Morse word can be defined as the number
        of bits in the 2-complement representation of the position
        modulo 2 which is what this function uses.  The running time
        is `O(\log n)` where `n` is the position desired.

        The `n`-th digit of the Generalized Thue Morse word can be defined as
        the sum of the digits of `n` written in the given base mod `m`,
        where `m` is the length of the given alphabet.

        INPUT:

        - ``n`` -- integer; the position
        - ``alphabet`` -- an alphabet (default: (0, 1)) of size at least 2
        - ``base`` -- integer (default: 2) greater than or equal to 2

        OUTPUT:

        0 or 1 -- the digit at the position
        letter -- the letter of alphabet at the position

        TESTS::

            sage: from sage.combinat.words.word_generators import WordGenerator
            sage: WordGenerator()._ThueMorseWord_nth_digit(0)
            0
            sage: WordGenerator()._ThueMorseWord_nth_digit(3)
            0
            sage: WordGenerator()._ThueMorseWord_nth_digit(32)
            1
            sage: WordGenerator()._ThueMorseWord_nth_digit(6, 'abc', base = 7)
            'a'

        Negative input::

            sage: words._ThueMorseWord_nth_digit(-7)
            Traceback (most recent call last):
            ...
            NotImplementedError: nth digit of Thue-Morse word is not implemented for negative value of n
        """
        if n < 0:
            raise NotImplementedError("nth digit of Thue-Morse word is not implemented for negative value of n")
        m = len(alphabet)
        if base == 2 and m == 2:
            for tn in count():
                if n == 0:
                    return alphabet[tn & 1]
                n &= n - 1
        elif base < 2 or m < 2:
            raise ValueError("base (=%s) and len(alphabet) (=%s) must be at least 2" % (base, m))
        else:
            return alphabet[ZZ(sum(ZZ(n).digits(base=base))).mod(m)]

    def FibonacciWord(self, alphabet=(0, 1), construction_method='recursive'):
        r"""
        Return the Fibonacci word on the given two-letter alphabet.

        INPUT:

        - ``alphabet`` -- any container of length two that is suitable to
          build an instance of OrderedAlphabet (``list``, ``tuple``, ``str``, ...)

        - ``construction_method`` -- can be any of the following:
          "recursive", "fixed point", "function" (see below for definitions)

        Recursive construction: the Fibonacci word is the limit of the
        following sequence of words: `S_0 = 0`, `S_1 = 01`,
        `S_n = S_{n-1} S_{n-2}` for `n \geq 2`.

        Fixed point construction: the Fibonacci word is the fixed point of the
        morphism: `0 \mapsto 01` and `1 \mapsto 0`. Hence, it can be constructed
        by the following read-write process:

        #. beginning at the first letter of `01`,
        #. if the next letter is `0`, append `01` to the word;
        #. if the next letter is `1`, append `1` to the word;
        #. move to the next letter of the word.

        Function: Over the alphabet `\{1, 2\}`, the `n`-th letter of the
        Fibonacci word is
        `\lfloor (n+2) \varphi \rfloor - \lfloor (n+1) \varphi \rfloor`
        where `\varphi=(1+\sqrt{5})/2` is the golden ratio.

        EXAMPLES::

            sage: w = words.FibonacciWord(construction_method='recursive'); w
            word: 0100101001001010010100100101001001010010...

        ::

            sage: v = words.FibonacciWord(construction_method='recursive', alphabet='ab'); v
            word: abaababaabaababaababaabaababaabaababaaba...

        ::

            sage: u = words.FibonacciWord(construction_method="fixed point"); u
            word: 0100101001001010010100100101001001010010...

        ::

            sage: words.FibonacciWord(construction_method="fixed point", alphabet=[4, 1])
            word: 4144141441441414414144144141441441414414...

        ::

            sage: words.FibonacciWord([0,1], 'function')                                # needs sage.symbolic
            word: 0100101001001010010100100101001001010010...
            sage: words.FibonacciWord('ab', 'function')                                 # needs sage.symbolic
            word: abaababaabaababaababaabaababaabaababaaba...

        TESTS::

            sage: from math import floor, sqrt
            sage: golden_ratio = (1 + sqrt(5))/2.0
            sage: a = golden_ratio / (1  + 2*golden_ratio)
            sage: wn = lambda n : int(floor(a*(n+2)) - floor(a*(n+1)))
            sage: f = Words([0,1])(wn); f
            word: 0100101001001010010100100101001001010010...
            sage: f[:10000] == w[:10000]
            True
            sage: f[:10000] == u[:10000]  #long time
            True
            sage: words.FibonacciWord("abc")
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements
        """
        W = InfiniteWords(alphabet)
        alphabet = W.alphabet()
        if alphabet.cardinality() != 2:
            raise TypeError("alphabet does not contain two distinct elements")
        a,b = alphabet

        if construction_method == "recursive":
            w = W(self._FibonacciWord_RecursiveConstructionIterator(alphabet),
                  datatype='iter')
            return w

        elif construction_method in ("fixed point", "fixed_point"):
            d = {b:[a],a:[a,b]}
            w = self.FixedPointOfMorphism(d, a)
            return w

        elif construction_method == "function":
            from sage.functions.other import floor
            from sage.misc.functional import sqrt
            phi = (1 + sqrt(5))/2 # the golden ratio
            f = lambda n:a if floor((n+2)*phi) - floor((n+1)*phi) == 2 else b
            return W(f)

        else:
            raise NotImplementedError

    def _FibonacciWord_RecursiveConstructionIterator(self,alphabet=(0,1)):
        r"""
        Iterate over the symbols of the Fibonacci word, as defined by
        the following recursive construction: the Fibonacci word is the
        limit of the sequence `S_0 = 0`, `S_1 = 01`, `S_n = S_{n-1}
        S_{n-2}` for `n \geq 2`.

        TESTS::

            sage: from sage.combinat.words.word_generators import WordGenerator
            sage: from itertools import islice
            sage: it = WordGenerator()._FibonacciWord_RecursiveConstructionIterator()
            sage: list(islice(it,13r))
            [0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1]
        """
        Fib0 = [0]
        Fib1 = [0,1]
        n = 0
        while True:
            it = iter(Fib1[n:])
            for i in it:
                n += 1
                yield alphabet[i]
            else:
                Fib1, Fib0 = Fib1 + Fib0, Fib1

    def FixedPointOfMorphism(self, morphism, first_letter):
        r"""
        Return the fixed point of the morphism beginning with
        ``first_letter``.

        A *fixed point* of a morphism `\varphi` is a word `w` such that
        `\varphi(w) = w`.

        INPUT:

        - ``morphism`` -- endomorphism prolongable on ``first_letter``. It
          must be something that WordMorphism's constructor understands
          (dict, str, ...).

        - ``first_letter`` -- the first letter of the fixed point

        OUTPUT: the fixed point of the morphism beginning with ``first_letter``

        EXAMPLES::

            sage: mu = {0:[0,1], 1:[1,0]}
            sage: tm = words.FixedPointOfMorphism(mu,0); tm
            word: 0110100110010110100101100110100110010110...
            sage: TM = words.ThueMorseWord()
            sage: tm[:1000] == TM[:1000]                                                # needs sage.modules
            True

        ::

            sage: mu = {0:[0,1], 1:[0]}
            sage: f = words.FixedPointOfMorphism(mu,0); f
            word: 0100101001001010010100100101001001010010...
            sage: F = words.FibonacciWord(); F
            word: 0100101001001010010100100101001001010010...
            sage: f[:1000] == F[:1000]                                                  # needs sage.modules
            True

        ::

            sage: fp = words.FixedPointOfMorphism('a->abc,b->,c->','a'); fp
            word: abc
        """
        return WordMorphism(morphism).fixed_point(letter=first_letter)

    def CodingOfRotationWord(self, alpha, beta, x=0, alphabet=(0,1)):
        r"""
        Return the infinite word obtained from the coding of rotation of
        parameters `(\alpha,\beta, x)` over the given two-letter alphabet.

        The *coding of rotation* corresponding to the parameters
        `(\alpha,\beta, x)` is the symbolic sequence `u = (u_n)_{n\geq 0}`
        defined over the binary alphabet `\{0, 1\}` by `u_n = 1` if
        `x+n\alpha\in[0, \beta[` and `u_n = 0` otherwise. See [AC03]_.

        EXAMPLES::

            sage: alpha = 0.45
            sage: beta = 0.48
            sage: words.CodingOfRotationWord(0.45, 0.48)
            word: 1101010101001010101011010101010010101010...

        ::

            sage: words.CodingOfRotationWord(0.45, 0.48, alphabet='xy')
            word: yyxyxyxyxyxxyxyxyxyxyyxyxyxyxyxxyxyxyxyx...

        TESTS::

            sage: words.CodingOfRotationWord(0.51,0.43,alphabet=[1,0,2])
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements
        """
        if len(set(alphabet)) != 2:
            raise TypeError("alphabet does not contain two distinct elements")
        from functools import partial
        f = partial(self._CodingOfRotationWord_function,alpha=alpha,beta=beta,x=x,alphabet=alphabet)
        w = InfiniteWords(alphabet)(f, datatype='callable')
        return w

    def _CodingOfRotationWord_function(self, n, alpha, beta, x=0, alphabet=(0,1)):
        r"""
        Internal function that returns the symbol in position `n` of the
        coding of rotation word corresponding to the parameters `\alpha`,
        `\beta`, and `x`.

        TESTS::

            sage: alpha, beta = 0.45, 0.48
            sage: words._CodingOfRotationWord_function(3, alpha, beta)
            1
            sage: words._CodingOfRotationWord_function(10, alpha, beta)
            0
            sage: words._CodingOfRotationWord_function(17, alpha, beta)
            0
        """
        hauteur = x + n * alpha
        fracH = hauteur.frac()
        if fracH < 0:
            fracH += 1
        if 0 <= fracH < beta:
            return alphabet[1]
        else:
            return alphabet[0]

    @rename_keyword(cf='slope')
    def CharacteristicSturmianWord(self, slope, alphabet=(0, 1), bits=None):
        r"""
        Return the characteristic Sturmian word (also called standard
        Sturmian word) of given slope.

        Over a binary alphabet `\{a,b\}`, the characteristic Sturmian
        word `c_\alpha` of irrational slope `\alpha` is the infinite word
        satisfying `s_{\alpha,0} = ac_\alpha` and `s'_{\alpha,0} = bc_\alpha`,
        where `s_{\alpha,0}` and `s'_{\alpha,0}` are respectively the lower
        and upper mechanical words with slope `\alpha` and intercept `0`.
        Equivalently, for irrational `\alpha`,
        `c_\alpha = s_{\alpha,\alpha} = s'_{\alpha,\alpha}`.

        Let `\alpha = [0, d_1 + 1, d_2, d_3, \ldots]` be the continued
        fraction expansion of `\alpha`. It has been shown that the
        characteristic Sturmian word of slope `\alpha` is also the limit of
        the sequence: `s_0 = b, s_1 = a, \ldots, s_{n+1} = s_n^{d_n} s_{n-1}`
        for `n > 0`.

        See Section 2.1 of [Loth02]_ for more details.

        INPUT:

        - ``slope`` -- the slope of the word. It can be one of the following:

          - real number in `]0, 1[`

          - iterable over the continued fraction expansion of a real
            number in `]0, 1[`

        - ``alphabet`` -- any container of length two that is suitable to
          build an instance of OrderedAlphabet (``list``, ``tuple``, ``str``,
          ...)

        - ``bits`` -- integer (optional and considered only if ``slope`` is
          a real number); the number of bits to consider when computing the
          continued fraction

        OUTPUT: word

        ALGORITHM:

        Let `[0, d_1 + 1, d_2, d_3, \ldots]` be the continued fraction
        expansion of `\alpha`. Then, the characteristic Sturmian word of
        slope `\alpha` is the limit of the sequence: `s_0 = b`, `s_1 = a`
        and `s_{n+1} = s_n^{d_n} s_{n-1}` for `n > 0`.

        EXAMPLES:

        From real slope::

            sage: words.CharacteristicSturmianWord(1/golden_ratio^2)                    # needs sage.symbolic
            word: 0100101001001010010100100101001001010010...
            sage: words.CharacteristicSturmianWord(4/5)                                 # needs sage.rings.real_mpfr
            word: 11110
            sage: words.CharacteristicSturmianWord(5/14)                                # needs sage.rings.real_mpfr
            word: 01001001001001
            sage: words.CharacteristicSturmianWord(pi - 3)                              # needs sage.symbolic
            word: 0000001000000100000010000001000000100000...

        From an iterator of the continued fraction expansion of a real::

            sage: def cf():
            ....:   yield 0
            ....:   yield 2
            ....:   while True: yield 1
            sage: F = words.CharacteristicSturmianWord(cf()); F                         # needs sage.rings.real_mpfr
            word: 0100101001001010010100100101001001010010...
            sage: Fib = words.FibonacciWord(); Fib
            word: 0100101001001010010100100101001001010010...
            sage: F[:10000] == Fib[:10000]                                              # needs sage.rings.real_mpfr
            True

        The alphabet may be specified::

            sage: words.CharacteristicSturmianWord(cf(), 'rs')                          # needs sage.rings.real_mpfr
            word: rsrrsrsrrsrrsrsrrsrsrrsrrsrsrrsrrsrsrrsr...

        The characteristic sturmian word of slope `(\sqrt{3}-1)/2`::

            sage: words.CharacteristicSturmianWord((sqrt(3)-1)/2)                       # needs sage.symbolic
            word: 0100100101001001001010010010010100100101...

        The same word defined from the continued fraction expansion of
        `(\sqrt{3}-1)/2`::

            sage: from itertools import cycle, chain
            sage: it = chain([0], cycle([2, 1]))
            sage: words.CharacteristicSturmianWord(it)
            word: 0100100101001001001010010010010100100101...

        The first terms of the standard sequence of the characteristic
        sturmian word of slope `(\sqrt{3}-1)/2`::

            sage: words.CharacteristicSturmianWord([0,2])
            word: 01
            sage: words.CharacteristicSturmianWord([0,2,1])
            word: 010
            sage: words.CharacteristicSturmianWord([0,2,1,2])
            word: 01001001
            sage: words.CharacteristicSturmianWord([0,2,1,2,1])
            word: 01001001010
            sage: words.CharacteristicSturmianWord([0,2,1,2,1,2])
            word: 010010010100100100101001001001
            sage: words.CharacteristicSturmianWord([0,2,1,2,1,2,1])
            word: 0100100101001001001010010010010100100101...

        TESTS::

            sage: words.CharacteristicSturmianWord([1,1,1],'xyz')
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain two distinct elements

        ::

            sage: words.CharacteristicSturmianWord(5/4)
            Traceback (most recent call last):
            ...
            ValueError: the argument slope (=5/4) must be in ]0,1[

        ::

            sage: words.CharacteristicSturmianWord(1/golden_ratio^2)                    # needs sage.symbolic
            word: 0100101001001010010100100101001001010010...
            sage: _.length()                                                            # needs sage.symbolic
            +Infinity

        ::

            sage: # needs sage.symbolic
            sage: a = words.LowerMechanicalWord(1/pi)[1:]
            sage: b = words.UpperMechanicalWord(1/pi)[1:]
            sage: c = words.CharacteristicSturmianWord(1/pi)
            sage: n = 500; a[:n] == b[:n] == c[:n]
            True

        ::

            sage: alpha = random()
            sage: c = words.CharacteristicSturmianWord(alpha)
            sage: l = words.LowerMechanicalWord(alpha)[1:]
            sage: u = words.UpperMechanicalWord(alpha)[1:]
            sage: i = 10000; j = i + 500; c[i:j] == l[i:j] == u[i:j]
            True

        ::

            sage: a, b = 207, 232
            sage: u = words.ChristoffelWord(a, b)
            sage: v = words.CharacteristicSturmianWord(a/(a+b))
            sage: v.length()
            439
            sage: u[1:-1] == v[:-2]
            True
        """
        if len(set(alphabet)) != 2:
            raise TypeError("alphabet does not contain two distinct elements")

        if slope in RR:
            if not 0 < slope < 1:
                msg = "the argument slope (=%s) must be in ]0,1[" % slope
                raise ValueError(msg)
            from sage.rings.continued_fraction import continued_fraction
            cf = continued_fraction(slope)
            if cf.length() == Infinity:
                parent = InfiniteWords(alphabet)
            else:
                parent = FiniteWords(alphabet)
            cf = iter(cf)
        elif isinstance(slope, Iterable):
            cf = iter(slope)
            parent = InfiniteWords(alphabet)
        else:
            raise TypeError("slope (=%s) must be a real number" % slope +
                            "or an iterable")
        w = parent(self._CharacteristicSturmianWord_LetterIterator(cf,alphabet),
                   datatype='iter')
        return w

    def _CharacteristicSturmianWord_LetterIterator(self, cf, alphabet=(0,1)):
        r"""
        Return an iterator over the symbols of the characteristic
        Sturmian word of slope ``cf``.

        INPUT:

        - ``cf`` -- iterator, the continued fraction expansion of a real
          number in `]0, 1[`
        - ``alphabet`` -- the alphabet (default: ``(0,1)``) of
          the output

        OUTPUT: iterator of letters

        ALGORITHM:

        Let `[0, d_1 + 1, d_2, d_3, \ldots]` be the continued fraction
        expansion of `\alpha`. Then, the characteristic Sturmian word of
        slope `\alpha` is the limit of the sequence: `s_0 = 1`, `s_1 = 0`
        and `s_{n+1} = s_n^{d_n} s_{n-1}` for `n > 0`.

        EXAMPLES::

            sage: continued_fraction(1/golden_ratio^2)[:8]                              # needs sage.symbolic
            [0; 2, 1, 1, 1, 1, 2]
            sage: cf = iter(_)                                                          # needs sage.symbolic
            sage: Word(words._CharacteristicSturmianWord_LetterIterator(cf))            # needs sage.symbolic
            word: 0100101001001010010100100101001010

        ::

            sage: # needs sage.symbolic
            sage: alpha = (sqrt(3)-1)/2
            sage: continued_fraction(alpha)[:10]
            [0; 2, 1, 2, 1, 2, 1, 2, 1, 2]
            sage: cf = iter(_)
            sage: Word(words._CharacteristicSturmianWord_LetterIterator(cf))
            word: 0100100101001001001010010010010100100101...
        """
        try:
            if next(cf) != 0:
                raise ValueError("the first term of the continued fraction expansion must be zero")
        except StopIteration:
            return

        s0 = [1]
        s1 = [0]
        try:
            e = next(cf)
        except StopIteration:
            return

        if not e >= 1:
            raise ValueError("the second term of the continued fraction expansion must be larger or equal to 1")
        s1, s0 = s1*(e-1) + s0, s1
        n = 0
        while True:
            try:
                for i in s1[n:]:
                    n += 1
                    yield alphabet[i]
                else:
                    s1, s0 = s1*next(cf) + s0, s1
            except StopIteration:
                return

    def KolakoskiWord(self, alphabet=(1,2)):
        r"""
        Return the Kolakoski word over the given alphabet and
        starting with the first letter of the alphabet.

        Let `A = \{a,b\}` be an alphabet, where `a` and `b` are two
        distinct positive integers. The Kolakoski word `K_{a,b}`
        over `A` and starting with `a` is the unique infinite word `w`
        such that `w = \Delta(w)`, where `\Delta(w)` is the word
        encoding the runs of `w` (see ``delta()`` method on words for
        more details).

        Note that `K_{a,b} \neq K_{b,a}`. On the other hand, the
        words `K_{a,b}` and `K_{b,a}` are the unique two words over `A`
        that are fixed by `\Delta`.

        Also note that the Kolakoski word is also known as the
        Oldenburger word.

        INPUT:

        - ``alphabet`` -- (default: (1,2)) an iterable of two positive
          integers

        OUTPUT: infinite word

        EXAMPLES:

        The usual Kolakoski word::

            sage: w = words.KolakoskiWord()
            sage: w
            word: 1221121221221121122121121221121121221221...
            sage: w.delta()
            word: 1221121221221121122121121221121121221221...

        The other Kolakoski word on the same alphabet::

            sage: w = words.KolakoskiWord(alphabet = (2,1))
            sage: w
            word: 2211212212211211221211212211211212212211...
            sage: w.delta()
            word: 2211212212211211221211212211211212212211...

        It is naturally generalized to any two integers alphabet::

            sage: w = words.KolakoskiWord(alphabet = (2,5))
            sage: w
            word: 2255222225555522552255225555522222555552...
            sage: w.delta()
            word: 2255222225555522552255225555522222555552...

        TESTS::

            sage: for i in range(1,10):
            ....:     for j in range(1,10):
            ....:         if i != j:
            ....:             w = words.KolakoskiWord(alphabet=(i,j))
            ....:             assert w[:50] == w.delta()[:50]

        ::

            sage: words.KolakoskiWord((0, 2))
            Traceback (most recent call last):
            ...
            ValueError: the alphabet (=(0, 2)) must consist of two distinct positive integers

        REFERENCES:

        .. [Kolakoski66] William Kolakoski, proposal 5304, American Mathematical Monthly
           72 (1965), 674; for a partial solution, see "Self Generating Runs,"
           by Necdet Üçoluk, Amer. Math. Mon. 73 (1966), 681-2.
        """
        a, b = alphabet
        if a not in ZZ or a <= 0 or b not in ZZ or b <= 0 or a == b:
            msg = 'the alphabet (=%s) must consist of two distinct positive integers' % (alphabet,)
            raise ValueError(msg)
        return InfiniteWords(alphabet)(self._KolakoskiWord_iterator(a, b), datatype='iter')

    def _KolakoskiWord_iterator(self, a=1, b=2):
        r"""
        Return an iterator over the Kolakoski word over ``{a,b}``
        and starting with ``a``.

        Let `A = \{a,b\}` be an alphabet, where `a` and `b` are two
        distinct positive integers. The Kolakoski word `K_{a,b}`
        over `A` and starting with `a` is the unique infinite word `w`
        such that `w = \Delta(w)`, where `\Delta(w)` is the word
        encoding the runs of `w` (see ``delta()`` method on words for
        more details).

        Note that `K_{a,b} \neq K_{b,a}`. On the other hand, the
        words `K_{a,b}` and `K_{b,a}` are the unique two words over `A`
        that are fixed by `\Delta`.

        INPUT:

        - ``a`` -- positive integer (default: 1); the first letter occurring
          in the returned Kolakoski word
        - ``b`` -- positive integer (default: 2); the second and last letter
          occurring in the returned Kolakoski word

        OUTPUT: iterator

        EXAMPLES:

        The first ten letters of `K_{3,5}`::

            sage: iter = words._KolakoskiWord_iterator(3, 5)
            sage: Word(iter)[:10]
            word: 3335553335

        See ``words.KolakoskiWord()`` for more documentation.
        """
        # First, we need to treat the basis case
        w = [a] * a
        for _ in range(a):
            yield a
        if a == 1:
            w.extend([b] * b)
            for _ in range(b):
                yield b
            w.pop(0)
        w.pop(0)
        # Letters swap function
        bar = lambda x: a if x == b else b
        current_letter = bar(w[-1])
        # Now we are ready to go in the recursive part
        while True:
            for _ in range(w[0]):
                yield current_letter
                w.append(current_letter)
            w.pop(0)
            current_letter = bar(current_letter)

    def LowerMechanicalWord(self, alpha, rho=0, alphabet=None):
        r"""
        Return the lower mechanical word with slope `\alpha` and
        intercept `\rho`

        The lower mechanical word `s_{\alpha,\rho}` with
        slope `\alpha` and intercept `\rho` is defined by
        `s_{\alpha,\rho}(n) = \lfloor\alpha(n+1) + \rho\rfloor -
        \lfloor\alpha n + \rho\rfloor`. [Loth02]_

        INPUT:

        - ``alpha`` -- real number such that `0 \leq\alpha\leq 1`

        - ``rho`` -- real number (default: 0)

        - ``alphabet`` -- iterable of two elements or ``None``
          (default: ``None``)

        OUTPUT: infinite word

        EXAMPLES::

            sage: words.LowerMechanicalWord(1/golden_ratio^2)                           # needs sage.symbolic
            word: 0010010100100101001010010010100100101001...
            sage: words.LowerMechanicalWord(1/5)                                        # needs sage.symbolic
            word: 0000100001000010000100001000010000100001...
            sage: words.LowerMechanicalWord(1/pi)                                       # needs sage.symbolic
            word: 0001001001001001001001000100100100100100...

        TESTS::

            sage: m = words.LowerMechanicalWord(1/golden_ratio^2)[1:]                   # needs sage.symbolic
            sage: s = words.CharacteristicSturmianWord(1/golden_ratio^2)                # needs sage.symbolic
            sage: m[:500] == s[:500]                                                    # needs sage.symbolic
            True

        Check that this returns a word in an alphabet (:issue:`10054`)::

            sage: words.UpperMechanicalWord(1/golden_ratio^2).parent()                  # needs sage.symbolic
            Infinite words over {0, 1}
        """
        if not 0 <= alpha <= 1:
            raise ValueError("parameter alpha (=%s) must be in [0,1]" % alpha)

        from sage.functions.other import floor
        from sage.combinat.words.alphabet import build_alphabet
        if alphabet is None or alphabet in ((0, 1), [0, 1]):
            alphabet = build_alphabet([0, 1])
            s = lambda n: floor(alpha*(n+1) + rho) - floor(alpha*n + rho)
        else:
            alphabet = build_alphabet(alphabet)
            card = alphabet.cardinality()
            if card != 2:
                raise TypeError("size of alphabet (=%s) must be two" % card)
            s = lambda n: alphabet[floor(alpha*(n+1) + rho) - floor(alpha*n + rho)]
        return InfiniteWords(alphabet)(s)

    def UpperMechanicalWord(self, alpha, rho=0, alphabet=None):
        r"""
        Return the upper mechanical word with slope `\alpha` and
        intercept `\rho`

        The upper mechanical word `s'_{\alpha,\rho}` with
        slope `\alpha` and intercept `\rho` is defined by
        `s'_{\alpha,\rho}(n) = \lceil\alpha(n+1) + \rho\rceil -
        \lceil\alpha n + \rho\rceil`. [Loth02]_

        INPUT:

        - ``alpha`` -- real number such that `0 \leq\alpha\leq 1`

        - ``rho`` -- real number (default: 0)

        - ``alphabet`` -- iterable of two elements or ``None``
          (default: ``None``)

        OUTPUT: infinite word

        EXAMPLES::

            sage: words.UpperMechanicalWord(1/golden_ratio^2)                           # needs sage.symbolic
            word: 1010010100100101001010010010100100101001...
            sage: words.UpperMechanicalWord(1/5)                                        # needs sage.symbolic
            word: 1000010000100001000010000100001000010000...
            sage: words.UpperMechanicalWord(1/pi)                                       # needs sage.symbolic
            word: 1001001001001001001001000100100100100100...

        TESTS::

            sage: m = words.UpperMechanicalWord(1/golden_ratio^2)[1:]                   # needs sage.symbolic
            sage: s = words.CharacteristicSturmianWord(1/golden_ratio^2)                # needs sage.symbolic
            sage: m[:500] == s[:500]                                                    # needs sage.symbolic
            True

        Check that this returns a word in an alphabet (:issue:`10054`)::

            sage: words.UpperMechanicalWord(1/golden_ratio^2).parent()                  # needs sage.symbolic
            Infinite words over {0, 1}
        """
        if not 0 <= alpha <= 1:
            raise ValueError("parameter alpha (=%s) must be in [0,1]" % alpha)

        from sage.functions.other import ceil
        from sage.combinat.words.alphabet import build_alphabet
        if alphabet is None or alphabet in ((0, 1), [0, 1]):
            alphabet = build_alphabet([0, 1])
            s = lambda n: ceil(alpha*(n+1) + rho) - ceil(alpha*n + rho)
        else:
            alphabet = build_alphabet(alphabet)
            card = alphabet.cardinality()
            if card != 2:
                raise TypeError("size of alphabet (=%s) must be two" % card)
            s = lambda n: alphabet[ceil(alpha*(n+1) + rho) - ceil(alpha*n + rho)]
        return InfiniteWords(alphabet)(s)

    def StandardEpisturmianWord(self, directive_word):
        r"""
        Return the standard episturmian word (or epistandard word) directed by
        directive_word. Over a 2-letter alphabet, this function
        gives characteristic Sturmian words.

        An infinite word `w` over a finite alphabet `A` is said to be
        *standard episturmian* (or *epistandard*) iff there exists an
        infinite word `x_1x_2x_3\cdots` over `A` (called the *directive
        word* of `w`) such that `w` is the limit as `n` goes to infinity of
        `Pal(x_1\cdots x_n)`, where `Pal` is the iterated palindromic closure
        function.

        Note that an infinite word is *episturmian* if it has the same set
        of factors as some epistandard word.

        See for instance [DJP2001]_, [JP2002]_, and [GJ2007]_.

        INPUT:

        - ``directive_word`` -- an infinite word or a period of a periodic
          infinite word

        EXAMPLES::

            sage: Fibonacci = words.StandardEpisturmianWord(Words('ab')('ab')); Fibonacci
            word: abaababaabaababaababaabaababaabaababaaba...
            sage: Tribonacci = words.StandardEpisturmianWord(Words('abc')('abc')); Tribonacci
            word: abacabaabacababacabaabacabacabaabacababa...
            sage: S = words.StandardEpisturmianWord(Words('abcd')('aabcabada')); S
            word: aabaacaabaaabaacaabaabaacaabaaabaacaabaa...
            sage: S = words.StandardEpisturmianWord(Fibonacci); S
            word: abaabaababaabaabaababaabaababaabaabaabab...
            sage: S[:25]
            word: abaabaababaabaabaababaaba
            sage: S = words.StandardEpisturmianWord(Tribonacci); S
            word: abaabacabaabaabacabaababaabacabaabaabaca...
            sage: words.StandardEpisturmianWord(123)
            Traceback (most recent call last):
            ...
            TypeError: directive_word is not a word, so it cannot be used to build an episturmian word
            sage: words.StandardEpisturmianWord(Words('ab'))
            Traceback (most recent call last):
            ...
            TypeError: directive_word is not a word, so it cannot be used to build an episturmian word
        """
        if not isinstance(directive_word, Word_class):
            raise TypeError("directive_word is not a word, so it cannot be used to build an episturmian word")
        epistandard = directive_word.parent()(
            self._StandardEpisturmianWord_LetterIterator(directive_word),
            datatype='iter')
        return epistandard

    def _StandardEpisturmianWord_LetterIterator(self, directive_word):
        r"""
        Internal iterating over the symbols of the standard episturmian
        word defined by the (directive) word directive_word.

        An infinite word `w` over a finite alphabet `A` is standard episturmian
        (or epistandard) iff there exists an infinite word `x_1x_2x_3\ldots`
        over `A` (called the directive word of `w`) such that `w` is the limit
        as `n` goes to infinity of `Pal(x_1x_2\cdots x_n)`, where `Pal` is the
        iterated palindromic closure function.

        INPUT:

        - ``directive_word`` -- an infinite word or a finite word; if
          directive_word is finite, then it is repeated to give
          an infinite word

        TESTS::

            sage: import itertools
            sage: it = words._StandardEpisturmianWord_LetterIterator(Word('ab'))
            sage: list(itertools.islice(it, 13r))
            ['a', 'b', 'a', 'a', 'b', 'a', 'b', 'a', 'a', 'b', 'a', 'a', 'b']
        """
        if isinstance(directive_word, FiniteWord_class):
            d = cycle(directive_word)
        else:
            d = iter(directive_word)
        W = directive_word.parent()
        w = W(next(d))
        n = 0
        while True:
            for x in w[n:]:
                n += 1
                yield x
            else:
                w = W(w * W(next(d))).palindromic_closure()

    def MinimalSmoothPrefix(self, n):
        r"""
        This function finds and returns the minimal smooth prefix of length
        ``n``.

        See [BMP2007]_ for a definition.

        INPUT:

        - ``n`` -- the desired length of the prefix

        OUTPUT: word; the prefix

        .. NOTE::

            Be patient, this function can take a really long time if asked
            for a large prefix.

        EXAMPLES::

            sage: words.MinimalSmoothPrefix(10)
            word: 1212212112
        """
        tab = []
        W = FiniteWords([1, 2])
        suff1 = W([1, 2, 2]).phi_inv()
        suff2 = W([2, 2]).phi_inv()
        w = [1]
        tab = _build_tab(1, tab, W)
        for k in range(1, n):
            if suff1._phi_inv_tab(tab) < suff2._phi_inv_tab(tab):
                w.append(1)
                tab = _build_tab(1, tab, W)
            else:
                w.append(2)
                tab = _build_tab(2, tab, W)
        return W(w)

    def RandomWord(self, n, m=2, alphabet=None):
        r"""
        Return a random word of length `n` over the given `m`-letter
        alphabet.

        INPUT:

        - ``n`` -- integer; the length of the word
        - ``m`` -- integer (default: 2); the size of the output alphabet
        - ``alphabet`` -- (default: `\{0,1,...,m-1\}`) any container of
          length m that is suitable to build an instance of
          OrderedAlphabet (``list``, ``tuple``, ``str``, ...)

        EXAMPLES::

            sage: words.RandomWord(10)         # random results
            word: 0110100101
            sage: words.RandomWord(10, 4)      # random results
            word: 0322313320
            sage: words.RandomWord(100, 7)     # random results
            word: 2630644023642516442650025611300034413310...
            sage: words.RandomWord(100, 7, range(-3,4))  # random results
            word: 1,3,-1,-1,3,2,2,0,1,-2,1,-1,-3,-2,2,0,3,0,-3,0,3,0,-2,-2,2,0,1,-3,2,-2,-2,2,0,2,1,-2,-3,-2,-1,0,...
            sage: words.RandomWord(100, 5, "abcde") # random results
            word: acebeaaccdbedbbbdeadeebbdeeebeaaacbadaac...
            sage: words.RandomWord(17, 5, "abcde")     # random results
            word: dcacbbecbddebaadd

        TESTS::

            sage: words.RandomWord(2,3,"abcd")
            Traceback (most recent call last):
            ...
            TypeError: alphabet does not contain 3 distinct elements
        """
        if alphabet is None:
            alphabet = list(range(m))
        if len(set(alphabet)) != m:
            raise TypeError("alphabet does not contain %s distinct elements" % m)
        return FiniteWords(alphabet)([alphabet[randint(0,m-1)] for i in range(n)])

    LowerChristoffelWord = LowerChristoffelWord

    ChristoffelWord = LowerChristoffelWord

    def UpperChristoffelWord(self, p, q, alphabet=(0,1)):
        r"""
        Return the upper Christoffel word of slope `p/q`, where
        `p` and `q` are relatively prime nonnegative
        integers, over the given alphabet.

        The *upper Christoffel word of slope `p/q`* is equal to the
        reversal of the lower Christoffel word of slope `p/q`.
        Equivalently, if `xuy` is the lower Christoffel word of
        slope `p/q`, where `x` and `y` are letters,
        then `yux` is the upper Christoffel word of slope
        `p/q` (because `u` is a palindrome).

        INPUT:

        - ``alphabet`` -- any container of length two that is suitable to build
          an instance of OrderedAlphabet (``list``, ``tuple``, ``str``, ...)

        EXAMPLES::

            sage: words.UpperChristoffelWord(1,0)
            word: 1

        ::

            sage: words.UpperChristoffelWord(0,1)
            word: 0

        ::

            sage: words.UpperChristoffelWord(1,1)
            word: 10

        ::

            sage: words.UpperChristoffelWord(4,7)
            word: 10100100100

        TESTS::

            sage: words.UpperChristoffelWord(51,43,"abc")
            Traceback (most recent call last):
            ...
            ValueError: alphabet must contain exactly two distinct elements
        """
        w = words.LowerChristoffelWord(p, q, alphabet=alphabet).reversal()
        return w

    @cached_method
    def _fibonacci_tile(self, n, q_0=None, q_1=3):
        r"""
        Return the word `q_n` defined by the recurrence below.

        The sequence `(q_n)_{n\in\NN}` is defined by `q_0=\varepsilon`,
        `q_1=3` and

        .. MATH::

            q_n = \begin{cases}
            q_{n-1}q_{n-2}       & \text{if} n\equiv 2 \mod 3, \\
            q_{n-1}\bar{q_{n-2}} & \text{if} n\equiv 0,1 \mod 3.
            \end{cases}

        where the operator `\bar{\,}` exchanges the `1` and `3`.

        INPUT:

        - ``n`` -- nonnegative integer
        - ``q_0`` -- first initial value (default: ``None``); it can be
          ``None``, 0, 1, 2 or 3
        - ``q_1`` -- second initial value (default: 3); it can be ``None``, 0,
          1, 2 or 3

        EXAMPLES::

            sage: for i in range(10): words._fibonacci_tile(i)
            word:
            word: 3
            word: 3
            word: 31
            word: 311
            word: 31131
            word: 31131133
            word: 3113113313313
            word: 311311331331331131133
            word: 3113113313313311311331331331131131

        REFERENCES:

        [BmBGL09]_
        """
        from sage.combinat.words.morphism import WordMorphism
        W = FiniteWords([0,1,2,3])
        bar = WordMorphism({0:0,1:3,3:1,2:2},codomain=W)
        if n == 0:
            a = [] if q_0 is None else [q_0]
            return W(a)
        elif n == 1:
            b = [] if q_1 is None else [q_1]
            return W(b)
        elif n % 3 == 2:
            u = self._fibonacci_tile(n-1,q_0,q_1)
            v = self._fibonacci_tile(n-2,q_0,q_1)
            return u * v
        else:
            u = self._fibonacci_tile(n-1,q_0,q_1)
            v = bar(self._fibonacci_tile(n-2,q_0,q_1))
            return u * v

    def fibonacci_tile(self, n):
        r"""
        Return the `n`-th Fibonacci Tile [BmBGL09]_.

        EXAMPLES::

            sage: for i in range(3): words.fibonacci_tile(i)                            # needs sage.modules
            Path: 3210
            Path: 323030101212
            Path: 3230301030323212323032321210121232121010...
        """
        w = self._fibonacci_tile(3*n+1)
        w = w**4
        from sage.combinat.words.paths import WordPaths
        P = WordPaths([0,1,2,3])
        l = list(w.partial_sums(start=3,mod=4))
        return P(l)[:-1]

    def dual_fibonacci_tile(self, n):
        r"""
        Return the `n`-th dual Fibonacci Tile [BmBGL09]_.

        EXAMPLES::

            sage: for i in range(4): words.dual_fibonacci_tile(i)                       # needs sage.modules
            Path: 3210
            Path: 32123032301030121012
            Path: 3212303230103230321232101232123032123210...
            Path: 3212303230103230321232101232123032123210...
        """
        w = self._fibonacci_tile(3*n+1,3,3)
        w = w**4
        from sage.combinat.words.paths import WordPaths
        P = WordPaths([0,1,2,3])
        l = list(w.partial_sums(start=3,mod=4))
        return P(l)[:-1]

    def _s_adic_iterator(self, sequence, letters):
        r"""
        Return the iterator over the `s`-adic infinite word obtained from a
        sequence of morphisms applied on letters where the hypothesis of
        nested prefixes is used.

        DEFINITION (from [Fogg]_):

        Let `w` be a infinite word over an alphabet `A = A_0`. A
        standard representation of `w` is obtained from a sequence of
        substitutions `\sigma_k : A_{k+1} \to A_k` and a sequence of letters
        `a_k \in A_k` such that:

        .. MATH::

            \lim_{k\to\infty} \sigma_0 \circ \sigma_1 \circ \cdots
            \sigma_k(a_k).

        Given a set of substitutions `S`, we say that the representation is
        `S`-adic standard if the substitutions are chosen in `S`.

        INPUT:

        - ``sequence`` -- an iterable sequence of morphisms. It may be finite
          or infinite
        - ``letters`` -- an iterable  sequence of letters. The image of the
          (i+1)-th letter under the (i+1)-th morphism must start with the i-th
          letter

        OUTPUT: iterator of letters

        EXAMPLES:

        Let's define three morphisms and compute the first nested successive
        prefixes of the `s`-adic word::

            sage: m1 = WordMorphism('e->gh,f->hg')
            sage: m2 = WordMorphism('c->ef,d->e')
            sage: m3 = WordMorphism('a->cd,b->dc')
            sage: Word(words._s_adic_iterator([m1],'e'))
            word: gh
            sage: Word(words._s_adic_iterator([m1,m2],'ec'))
            word: ghhg
            sage: Word(words._s_adic_iterator([m1,m2,m3],'eca'))
            word: ghhggh

        If the letters don't satisfy the hypothesis of the algorithm, an
        error is raised::

            sage: Word(words._s_adic_iterator([m1,m2,m3],'ecb'))
            Traceback (most recent call last):
            ...
            ValueError: the hypothesis of the algorithm used is not satisfied; the image of the 3-th letter (=b) under the 3-th morphism (=a->cd, b->dc) should start with the 2-th letter (=c)

        Two examples of infinite `s`-adic words::

            sage: tm = WordMorphism('a->ab,b->ba')
            sage: fib = WordMorphism('a->ab,b->a')
            sage: from itertools import repeat
            sage: Word(words._s_adic_iterator(repeat(tm),repeat('a')))
            word: abbabaabbaababbabaababbaabbabaabbaababba...
            sage: Word(words._s_adic_iterator(repeat(fib),repeat('a')))
            word: abaababaabaababaababaabaababaabaababaaba...

        A less trivial infinite `s`-adic word::

            sage: D = {4:tm,5:fib}
            sage: tmword = words.ThueMorseWord([4,5])
            sage: it = (D[a] for a in tmword)
            sage: Word(words._s_adic_iterator(it, repeat('a')))
            word: abbaababbaabbaabbaababbaababbaabbaababba...

        The morphism `\sigma: a \mapsto ba, b \mapsto b` cannot satisfy the
        hypothesis of the algorithm (nested prefixes)::

            sage: sigma = WordMorphism('a->ba,b->b')
            sage: Word(words._s_adic_iterator(repeat(sigma),repeat('a')))
            Traceback (most recent call last):
            ...
            ValueError: the hypothesis of the algorithm used is not satisfied; the image of the 2-th letter (=a) under the 2-th morphism (=a->ba, b->b) should start with the 1-th letter (=a)

        AUTHORS:

        - Sébastien Labbé (2009-12-18): initial version
        """
        from itertools import tee
        sequence_it,sequence = tee(sequence)
        m = next(sequence_it)
        codomain = m.codomain()
        p = codomain.identity_morphism()
        letters_it,letters = tee(letters)
        precedent_letter = m(next(letters_it))[0]

        yield precedent_letter
        for (i,(m,a)) in enumerate(zip(sequence, letters)):
            if not precedent_letter == m(a)[0]:
                raise ValueError("the hypothesis of the algorithm used is not satisfied; the image of the %s-th letter (=%s) under the %s-th morphism (=%s) should start with the %s-th letter (=%s)" % (i+1,a,i+1,m,i,precedent_letter))
            w = p(m(a)[1:])
            yield from w
            p = p * m
            precedent_letter = a

    def s_adic(self, sequence, letters, morphisms=None):
        r"""
        Return the `s`-adic infinite word obtained from a sequence of
        morphisms applied on a letter.

        DEFINITION (from [Fogg]_):

        Let `w` be a infinite word over an alphabet `A = A_0`. A
        standard representation of `w` is obtained from a sequence of
        substitutions `\sigma_k : A_{k+1} \to A_k` and a sequence of letters
        `a_k \in A_k` such that:

        .. MATH::

            \lim_{k\to\infty} \sigma_0 \circ \sigma_1 \circ \cdots
            \sigma_k(a_k).

        Given a set of substitutions `S`, we say that the representation is
        `S`-adic standard if the substitutions are chosen in `S`.

        INPUT:

        - ``sequence`` -- an iterable sequence of indices or of morphisms. It
          may be finite or infinite. If ``sequence`` is infinite, the image
          of the `(i+1)`-th letter under the `(i+1)`-th morphism must start
          with the `i`-th letter.

        - ``letters`` -- a letter or a sequence of letters

        - ``morphisms`` -- dict, list, callable or ``None`` (default:
          ``None``) an object that maps indices to morphisms. If ``None``, then
          ``sequence`` must consist of morphisms.

        OUTPUT: a word

        EXAMPLES:

        Let us define three morphisms and compute the first nested successive
        prefixes of the `s`-adic word::

            sage: m1 = WordMorphism('e->gh,f->hg')
            sage: m2 = WordMorphism('c->ef,d->e')
            sage: m3 = WordMorphism('a->cd,b->dc')
            sage: words.s_adic([m1],'e')
            word: gh
            sage: words.s_adic([m1,m2],'ec')
            word: ghhg
            sage: words.s_adic([m1,m2,m3],'eca')
            word: ghhggh

        When the given sequence of morphism is finite, one may simply give
        the last letter, i.e. ``'a'``, instead of giving all of them,
        i.e. ``'eca'``::

            sage: words.s_adic([m1,m2,m3],'a')
            word: ghhggh
            sage: words.s_adic([m1,m2,m3],'b')
            word: ghghhg

        If the letters don't satisfy the hypothesis of the algorithm
        (nested prefixes), an error is raised::

            sage: words.s_adic([m1,m2,m3],'ecb')
            Traceback (most recent call last):
            ...
            ValueError: the hypothesis of the algorithm used is not satisfied; the image of the 3-th letter (=b) under the 3-th morphism (=a->cd, b->dc) should start with the 2-th letter (=c)

        Let's define the Thue-Morse morphism and the Fibonacci morphism
        which will be used below to illustrate more examples and let's import
        the ``repeat`` tool from the ``itertools``::

            sage: tm = WordMorphism('a->ab,b->ba')
            sage: fib = WordMorphism('a->ab,b->a')
            sage: from itertools import repeat

        Two trivial examples of infinite `s`-adic words::

            sage: words.s_adic(repeat(tm),repeat('a'))
            word: abbabaabbaababbabaababbaabbabaabbaababba...

        ::

            sage: words.s_adic(repeat(fib),repeat('a'))
            word: abaababaabaababaababaabaababaabaababaaba...

        A less trivial infinite `s`-adic word::

            sage: D = {4:tm,5:fib}
            sage: tmword = words.ThueMorseWord([4,5])
            sage: it = (D[a] for a in tmword)
            sage: words.s_adic(it, repeat('a'))
            word: abbaababbaabbaabbaababbaababbaabbaababba...

        The same thing using a sequence of indices::

            sage: tmword = words.ThueMorseWord([0,1])
            sage: words.s_adic(tmword, repeat('a'), [tm,fib])
            word: abbaababbaabbaabbaababbaababbaabbaababba...

        The correspondence of the indices may be given as a dict::

            sage: words.s_adic(tmword, repeat('a'), {0:tm,1:fib})
            word: abbaababbaabbaabbaababbaababbaabbaababba...

        because dict are more versatile for indices::

            sage: tmwordTF = words.ThueMorseWord('TF')
            sage: words.s_adic(tmwordTF, repeat('a'), {'T':tm,'F':fib})
            word: abbaababbaabbaabbaababbaababbaabbaababba...

        or by a callable::

            sage: f = lambda n: tm if n == 0 else fib
            sage: words.s_adic(words.ThueMorseWord(), repeat('a'), f)
            word: abbaababbaabbaabbaababbaababbaabbaababba...

        Random infinite `s`-adic words::

            sage: from sage.misc.prandom import randint
            sage: def it():
            ....:   while True: yield randint(0,1)
            sage: words.s_adic(it(), repeat('a'), [tm,fib])  # random
            word: abbaabababbaababbaabbaababbaabababbaabba...
            sage: words.s_adic(it(), repeat('a'), [tm,fib])  # random
            word: abbaababbaabbaababbaababbaabbaababbaabba...
            sage: words.s_adic(it(), repeat('a'), [tm,fib])  # random
            word: abaaababaabaabaaababaabaaababaaababaabaa...

        An example where the sequences cycle on two morphisms and two
        letters::

            sage: G = WordMorphism('a->cd,b->dc')
            sage: H = WordMorphism('c->ab,d->ba')
            sage: from itertools import cycle
            sage: words.s_adic([G,H],'ac')
            word: cddc
            sage: words.s_adic(cycle([G,H]),cycle('ac'))
            word: cddcdccddccdcddcdccdcddccddcdccddccdcddc...

        The morphism `\sigma: a\mapsto ba, b\mapsto b` can't satisfy the
        hypothesis of the nested prefixes, but one may compute arbitrarily
        long finite words having the limit `\sigma^\omega(a)`::

            sage: sigma = WordMorphism('a->ba,b->b')
            sage: words.s_adic(repeat(sigma),repeat('a'))
            Traceback (most recent call last):
            ...
            ValueError: the hypothesis of the algorithm used is not satisfied; the image of the 2-th letter (=a) under the 2-th morphism (=a->ba, b->b) should start with the 1-th letter (=a)
            sage: words.s_adic([sigma],'a')
            word: ba
            sage: words.s_adic([sigma,sigma],'a')
            word: bba
            sage: words.s_adic([sigma]*3,'a')
            word: bbba
            sage: words.s_adic([sigma]*4,'a')
            word: bbbba
            sage: words.s_adic([sigma]*5,'a')
            word: bbbbba
            sage: words.s_adic([sigma]*6,'a')
            word: bbbbbba
            sage: words.s_adic([sigma]*7,'a')
            word: bbbbbbba

        The following examples illustrates an `S`-adic word defined over an
        infinite set `S` of morphisms `x_h`::

            sage: x = lambda h:WordMorphism({1:[2],2:[3]+[1]*(h+1),3:[3]+[1]*h})
            sage: for h in [0,1,2,3]:
            ....:     print("{} {}".format(h, x(h)))
            0 1->2, 2->31, 3->3
            1 1->2, 2->311, 3->31
            2 1->2, 2->3111, 3->311
            3 1->2, 2->31111, 3->3111
            sage: w = Word(lambda n : valuation(n+1, 2) ); w
            word: 0102010301020104010201030102010501020103...
            sage: s = words.s_adic(w, repeat(3), x); s
            word: 3232232232322322322323223223232232232232...
            sage: prefixe = s[:10000]
            sage: list(map(prefixe.number_of_factors, range(15)))
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            sage: [_[i+1] - _[i] for i in range(len(_)-1)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        TESTS::

            sage: tm = WordMorphism('a->ab,b->ba')
            sage: fib = WordMorphism('a->ab,b->a')
            sage: w = words.s_adic([fib,tm,tm,fib,tm,fib]*3,'a')
            sage: w
            word: abaaabaababaabaaababaaababaaabaababaabaa...
            sage: w.length()
            32400
            sage: w.parent()
            Finite words over {'a', 'b'}
            sage: type(w)
            <class 'sage.combinat.words.word.FiniteWord_callable_with_caching'>

        ::

            sage: words.s_adic([fib,tm,tm,fib,tm,fib],'aaaaaaa')
            word: abaaabaababaabaaababaaababaaabaababa

        ::

            sage: words.s_adic([0,1,0,1,0,1,0,1],'a',[tm,fib])
            word: abbaabababbaabbaababbaababbaabababbaabba...

        ::

            sage: words.s_adic([fib,fib],'bb')
            Traceback (most recent call last):
            ...
            ValueError: the hypothesis of the algorithm used is not satisfied; the image of the 2-th letter (=b) under the 2-th morphism (=a->ab, b->a) should start with the 1-th letter (=b)

        Test on different letters::

            sage: tm = WordMorphism({0:[0,1], 1:[1,0]})
            sage: fib = WordMorphism({0:[0,1], 1:[0]})
            sage: f = lambda n: tm if n == 0 else fib
            sage: words.s_adic(words.ThueMorseWord(), repeat(0), f)
            word: 0110010110011001100101100101100110010110...

        Testing the message error for the third argument::

            sage: words.s_adic(words.ThueMorseWord(), repeat(0), 5)
            Traceback (most recent call last):
            ...
            TypeError: morphisms (=5) must be None, callable or provide a __getitem__ method

        AUTHORS:

        - Sébastien Labbé (2009-12-18): initial version
        """
        if morphisms is None:
            seq = sequence
        elif hasattr(morphisms, '__getitem__'):
            seq = (morphisms[i] for i in sequence)
        elif callable(morphisms):
            seq = (morphisms(i) for i in sequence)
        else:
            raise TypeError("morphisms (=%s) must be None, callable or provide a __getitem__ method" % morphisms)

        from sage.combinat.words.word import FiniteWord_class
        if isinstance(sequence,(tuple,list,str,FiniteWord_class)) \
        and hasattr(letters, "__len__") and len(letters) == 1:
            from sage.misc.misc_c import prod
            return prod(seq)(letters)

        from itertools import tee
        seq_it, seq = tee(seq)
        m = next(seq_it)
        W = m.codomain()

        kwds = {}
        kwds['data'] = self._s_adic_iterator(seq, letters)
        kwds['datatype'] = 'iter'
        kwds['caching'] = True
        #kwds['check'] = False
        return W.shift()(**kwds)

    def PalindromicDefectWord(self, k=1, alphabet='ab'):
        r"""
        Return the finite word `w = a b^k a b^{k-1} a a b^{k-1} a b^{k} a`.

        As described by Brlek, Hamel, Nivat and Reutenauer in [BHNR2004]_, this
        finite word `w` is such that the infinite periodic word `w^{\omega}`
        has palindromic defect ``k``.

        INPUT:

        - ``k`` -- positive integer (default: 1)

        - ``alphabet`` -- iterable of size two (default: ``'ab'``)

        OUTPUT: finite word

        EXAMPLES::

            sage: words.PalindromicDefectWord(10)
            word: abbbbbbbbbbabbbbbbbbbaabbbbbbbbbabbbbbbb...

        ::

            sage: w = words.PalindromicDefectWord(3)
            sage: w
            word: abbbabbaabbabbba
            sage: w.defect()
            0
            sage: (w^2).defect()
            3
            sage: (w^3).defect()
            3

        On other alphabets::

            sage: words.PalindromicDefectWord(3, alphabet='cd')
            word: cdddcddccddcdddc
            sage: words.PalindromicDefectWord(3, alphabet=['c', 3])
            word: c333c33cc33c333c

        TESTS::

            sage: k = 25
            sage: (words.PalindromicDefectWord(k)^2).defect()
            25

        If k is negative or zero, then we get the same word::

            sage: words.PalindromicDefectWord(0)
            word: aaaaaa
            sage: words.PalindromicDefectWord(-3)
            word: aaaaaa
        """
        kk = k-1
        a, b = alphabet
        if not (isinstance(a, str) and isinstance(b, str)):
            a, b = (a,), (b,)
        w = a + b*k + a + b*kk + a + a + b*kk + a + b*k + a
        return FiniteWords(alphabet)(w)

    def BaumSweetWord(self):
        r"""
        Return the Baum-Sweet Word.

        The Baum-Sweet Sequence is an infinite word over the alphabet `\{0,1\}`
        defined by the following string substitution rules:

        `00 \rightarrow 0000`

        `01 \rightarrow 1001`

        `10 \rightarrow 0100`

        `11 \rightarrow 1101`

        The substitution rule above can be considered as a morphism on the
        submonoid of `\{0,1\}` generated by `\{00,01,10,11\}` (which is a free
        monoid on these generators).

        It is also defined as the concatenation of the terms from the Baum-Sweet
        Sequence:

        .. MATH::

            b_n = \begin{cases}
            0,      & \text{if } n = 0 \\
            1,      & \text{if } m \text{ is even} \\
            b_{\frac{m-1}{2}}, & \text{if } m \text{ is odd}
            \end{cases}

        where `n=m4^k` and `m` is not divisible by 4 if `m \neq 0`.

        The individual terms of the Baum-Sweet Sequence are also given by:

        .. MATH::

            b_n = \begin{cases}
            1,      & \text{if the binary representation of} n \text{ contains no block of consecutive 0s of odd length}\\
            0,      & \text{otherwise}\\
            \end{cases}\\

        for `n > 0` with `b_0 = 1`.

        For more information see:
        :wikipedia:`Baum-Sweet_sequence`.

        EXAMPLES:

        Baum-Sweet Word::

            sage: w = words.BaumSweetWord(); w
            word: 1101100101001001100100000100100101001001...

        Block Definition::

            sage: w = words.BaumSweetWord()
            sage: f = lambda n: '1' if all(len(x)%2==0 for x in bin(n)[2:].split('1')) else '0'
            sage: all(f(i) == w[i] for i in range(1,100))
            True
        """
        outer = WordMorphism('a->00,b->01,c->10,d->11')
        inner = WordMorphism('a->aa,b->cb,c->ba,d->db')
        return outer(inner.fixed_point('d'))


words = WordGenerator()
