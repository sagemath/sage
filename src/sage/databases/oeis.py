# -*- coding: utf-8 -*-
r"""
The On-Line Encyclopedia of Integer Sequences (OEIS)

You can query the OEIS (Online Database of Integer Sequences) through Sage in
order to:

    - identify a sequence from its first terms.
    - obtain more terms, formulae, references, etc. for a given sequence.


AUTHORS:

- Thierry Monteil (2012-02-10 -- 2013-06-21): initial version.

- Vincent Delecroix (2014): modifies continued fractions because of :trac:`14567`

- Moritz Firsching (2016): modifies handling of dead sequence, see :trac:`17330`

- Thierry Monteil (2019): refactorization (unique representation :trac:`28480`,
  laziness :trac:`28627`)

EXAMPLES::

        sage: oeis
        The On-Line Encyclopedia of Integer Sequences (https://oeis.org/)

What about a sequence starting with `3, 7, 15, 1` ?

::

    sage: search = oeis([3, 7, 15, 1], max_results=4) ; search  # optional -- internet
    0: A001203: Simple continued fraction expansion of Pi.
    1: A240698: Partial sums of divisors of n, cf. A027750.
    2: A082495: a(n) = (2^n - 1) mod n.
    3: A165416: Irregular array read by rows: The n-th row contains those distinct positive integers that each, when written in binary, occurs as a substring in binary n.

    sage: [u.id() for u in search]                      # optional -- internet
    ['A001203', 'A240698', 'A082495', 'A165416']
    sage: c = search[0] ; c                             # optional -- internet
    A001203: Simple continued fraction expansion of Pi.

::

    sage: c.first_terms(15)                             # optional -- internet
    (3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1)

    sage: c.examples()                                  # optional -- internet
    0: Pi = 3.1415926535897932384...
    1:    = 3 + 1/(7 + 1/(15 + 1/(1 + 1/(292 + ...))))
    2:    = [a_0; a_1, a_2, a_3, ...] = [3; 7, 15, 1, 292, ...].

    sage: c.comments()                                  # optional -- internet
    0: The first 5821569425 terms were computed by _Eric W. Weisstein_ on Sep 18 2011.
    1: The first 10672905501 terms were computed by _Eric W. Weisstein_ on Jul 17 2013.
    2: The first 15000000000 terms were computed by _Eric W. Weisstein_ on Jul 27 2013.
    3: The first 30113021586 terms were computed by _Syed Fahad_ on Apr 27 2021.

::

    sage: x = c.natural_object() ; type(x)              # optional -- internet
    <class 'sage.rings.continued_fraction.ContinuedFraction_periodic'>

    sage: x.convergents()[:7]                           # optional -- internet
    [3, 22/7, 333/106, 355/113, 103993/33102, 104348/33215, 208341/66317]

    sage: RR(x.value())                                 # optional -- internet
    3.14159265358979
    sage: RR(x.value()) == RR(pi)                       # optional -- internet
    True

What about posets ? Are they hard to count ? To which other structures are they
related ?

::

    sage: [Posets(i).cardinality() for i in range(10)]
    [1, 1, 2, 5, 16, 63, 318, 2045, 16999, 183231]
    sage: oeis(_)                                       # optional -- internet
    0: A000112: Number of partially ordered sets ("posets") with n unlabeled elements.
    sage: p = _[0]                                      # optional -- internet

::

    sage: 'hard' in p.keywords()                        # optional -- internet
    True
    sage: len(p.formulas())                             # optional -- internet
    0
    sage: len(p.first_terms())                          # optional -- internet
    17

::

    sage: p.cross_references(fetch=True)        # optional -- internet   # random
    0: A000798: Number of different quasi-orders (or topologies, or transitive digraphs) with n labeled elements.
    1: A001035: Number of partially ordered sets ("posets") with n labeled elements (or labeled acyclic transitive digraphs).
    2: A001930: Number of topologies, or transitive digraphs with n unlabeled nodes.
    3: A006057: Number of topologies on n labeled points satisfying axioms T_0-T_4.
    4: A079263: Number of constrained mixed models with n factors.
    5: A079265: Number of antisymmetric transitive binary relations on n unlabeled points.
    6: A263859: Triangle read by rows: T(n,k) (n>=1, k>=0) is the number of posets with n elements and rank k (or depth k+1).
    7: A316978: Number of factorizations of n into factors > 1 with no equivalent primes.
    8: A319559: Number of non-isomorphic T_0 set systems of weight n.
    9: A326939: Number of T_0 sets of subsets of {1..n} that cover all n vertices.
    10: A326943: Number of T_0 sets of subsets of {1..n} that cover all n vertices and are closed under intersection.
    ...


What does the Taylor expansion of the `e^{e^x-1}` function have to do with
primes ?

::

    sage: x = var('x') ; f(x) = e^(e^x - 1)
    sage: L = [a*factorial(b) for a,b in taylor(f(x), x, 0, 20).coefficients()] ; L
    [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597,
    27644437, 190899322, 1382958545, 10480142147, 82864869804, 682076806159,
    5832742205057, 51724158235372]

    sage: oeis(L)                                       # optional -- internet
    0: A000110: Bell or exponential numbers: number of ways to partition a set of n labeled elements.
    1: A292935: E.g.f.: exp(exp(-x) - 1).

    sage: b = _[0]                                      # optional -- internet

    sage: b.formulas()[0]                               # optional -- internet
    'E.g.f.: exp(exp(x) - 1).'

    sage: [i for i in b.comments() if 'prime' in i][-1]     # optional -- internet
    'Number n is prime if mod(a(n)-2,n) = 0. -_Dmitry Kruchinin_, Feb 14 2012'

    sage: [n for n in range(2, 20) if (b(n)-2) % n == 0]    # optional -- internet
    [2, 3, 5, 7, 11, 13, 17, 19]


.. SEEALSO::

    - If you plan to do a lot of automatic searches for subsequences, you
      should consider installing :mod:`SloaneEncyclopedia
      <sage.databases.sloane>`, a local partial copy of the OEIS.
    - Some infinite OEIS sequences are implemented in Sage, via the
      :mod:`sloane_functions <sage.combinat.sloane_functions>` module.

.. TODO::

    - in case of flood, suggest the user to install the off-line database instead.
    - interface with the off-line database (or reimplement it).

Classes and methods
-------------------
"""

# ****************************************************************************
#       Copyright (C) 2012 Thierry Monteil <sage!lma.metelu.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from urllib.request import urlopen
from urllib.parse import urlencode
from ssl import SSLContext

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.cpython.string import bytes_to_str
from sage.rings.integer import Integer
from sage.misc.verbose import verbose
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.temporary_file import tmp_filename
from sage.misc.unknown import Unknown
from sage.misc.html import HtmlFragment
from sage.repl.preparse import preparse

from collections import defaultdict
import re

oeis_url = 'https://oeis.org/'


def _fetch(url):
    r"""
    Fetch the given ``url``.

    INPUT:

    - ``url`` -- a string corresponding to the URL to be fetched.

    OUTPUT:

    - a string representing the fetched web page.

    TESTS::

        sage: from sage.databases.oeis import _fetch, oeis_url
        sage: _fetch(oeis_url + 'hints.html')[-8:-1]            # optional -- internet
        '</html>'
    """
    try:
        verbose("Fetching URL %s ..." % url, caller_name='OEIS')
        f = urlopen(url, context=SSLContext())
        result = f.read()
        f.close()
        return bytes_to_str(result)
    except IOError as msg:
        raise IOError("%s\nError fetching %s." % (msg, url))


def _urls(html_string):
    r"""
    Return the list of URLs contained in ``html_string``.

    Only URLs provided by HTML hyperlinks (``href`` attribute of ``<a>`` tags)
    in are returned, not text strings starting with ``http://``.

    INPUT:

    - ``html_string`` -- a string representing some HTML code.

    OUTPUT:

    - a list of (string) URLs contained in ``html_string``.

    EXAMPLES::

        sage: from sage.databases.oeis import _urls
        sage: html = 'http://example.com is not a link, but <a href="http://sagemath.org/">sagemath</a> is'
        sage: _urls(html)
        ['http://sagemath.org/']

    """
    urls = []
    from html.parser import HTMLParser

    class MyHTMLParser(HTMLParser):
        def handle_starttag(self, tag, attrs):
            if tag == 'a':
                for attr in attrs:
                    if attr[0] == 'href':
                        urls.append(attr[1])
    MyHTMLParser().feed(html_string)
    return urls


def to_tuple(string):
    """
    Convert a string to a tuple of integers.

    EXAMPLES::

        sage: from sage.databases.oeis import to_tuple
        sage: to_tuple('12,55,273')
        (12, 55, 273)
    """
    return tuple(Integer(x) for x in string.split(",") if x)


class OEIS:
    r"""
    The On-Line Encyclopedia of Integer Sequences.

    ``OEIS`` is a class representing the On-Line Encyclopedia of Integer
    Sequences. You can query it using its methods, but ``OEIS`` can also be
    called directly with three arguments:

    - ``query`` - it can be:

      - a string representing an OEIS ID (e.g. 'A000045').
      - an integer representing an OEIS ID (e.g. 45).
      - a list representing a sequence of integers.
      - a string, representing a text search.

    - ``max_results`` - (integer, default: 30) the maximum number of
      results to return, they are sorted according to their relevance. In
      any cases, the OEIS website will never provide more than 100 results.

    - ``first_result`` - (integer, default: 0) allow to skip the
      ``first_result`` first results in the search, to go further.
      This is useful if you are looking for a sequence that may appear
      after the 100 first found sequences.

    OUTPUT:

    - if ``query`` is an integer or an OEIS ID (e.g. 'A000045'), returns
      the associated OEIS sequence.

    - if ``query`` is a string, returns a tuple of OEIS sequences whose
      description corresponds to the query. Those sequences can be used
      without the need to fetch the database again.

    - if ``query`` is a list or tuple of integers, returns a tuple of
      OEIS sequences containing it as a subsequence. Those sequences
      can be used without the need to fetch the database again.

    EXAMPLES::

        sage: oeis
        The On-Line Encyclopedia of Integer Sequences (https://oeis.org/)

    A particular sequence can be called by its A-number or number::

        sage: oeis('A000040')                           # optional -- internet
        A000040: The prime numbers.

        sage: oeis(45)                                  # optional -- internet
        A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

    The database can be searched by subsequence::

        sage: search = oeis([1,2,3,5,8,13]) ; search    # optional -- internet
        0: A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
        1: A290689: Number of transitive rooted trees with n nodes.
        2: A027926: Triangular array T read by rows: T(n,0) = T(n,2n) = 1 for n >= 0; T(n,1) = 1 for n >= 1; T(n,k) = T(n-1,k-2) + T(n-1,k-1) for k = 2..2n-1, n >= 2.

        sage: fibo = search[0]                          # optional -- internet

        sage: fibo.name()                               # optional -- internet
        'Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.'

        sage: print(fibo.first_terms())                 # optional -- internet
        (0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597,
        2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418,
        317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465,
        14930352, 24157817, 39088169, 63245986, 102334155)

        sage: fibo.cross_references()[0]                # optional -- internet
        'A039834'

        sage: fibo == oeis(45)                          # optional -- internet
        True

        sage: sfibo = oeis('A039834')
        sage: sfibo.first_terms()                       # optional -- internet
        (1, 1, 0, 1, -1, 2, -3, 5, -8, 13, -21, 34, -55, 89, -144, 233,
        -377, 610, -987, 1597, -2584, 4181, -6765, 10946, -17711, 28657,
        -46368, 75025, -121393, 196418, -317811, 514229, -832040, 1346269,
        -2178309, 3524578, -5702887, 9227465, -14930352, 24157817)

        sage: tuple(abs(i) for i in sfibo.first_terms())[2:20] == fibo.first_terms()[:18]   # optional -- internet
        True

        sage: fibo.formulas()[4]                        # optional -- internet
        'F(n) = F(n-1) + F(n-2) = -(-1)^n F(-n).'

        sage: fibo.comments()[6]                        # optional -- internet
        "F(n+2) = number of binary sequences of length n that have no
        consecutive 0's."

        sage: fibo.links()[0]                           # optional -- internet
        'https://oeis.org/A000045/b000045.txt'

    The database can be searched by description::

        sage: oeis('prime gap factorization', max_results=4) # optional --internet
        0: A073491: Numbers having no prime gaps in their factorization.
        1: A073485: Product of any number of consecutive primes; squarefree numbers with no gaps in their prime factorization.
        2: A073490: Number of prime gaps in factorization of n.
        3: A073492: Numbers having at least one prime gap in their factorization.

    .. NOTE::

        The following will fetch the OEIS database only once::

            sage: oeis([1,2,3,5,8,13])                  # optional -- internet
            0: A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
            1: A290689: Number of transitive rooted trees with n nodes.
            2: A027926: Triangular array T read by rows: T(n,0) = T(n,2n) = 1 for n >= 0; T(n,1) = 1 for n >= 1; T(n,k) = T(n-1,k-2) + T(n-1,k-1) for k = 2..2n-1, n >= 2.

            sage: oeis('A000045')                       # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

        Indeed, due to some caching mechanism, the sequence is not re-created
        when called from its ID.

    TESTS::

        sage: oeis((1,2,5,16,61))    # optional -- internet
        0: A000111: ...
    """

    def __call__(self, query, max_results=3, first_result=0):
        r"""
        See the documentation of :class:`OEIS`.

        TESTS::

            sage: oeis()
            Traceback (most recent call last):
            ...
            TypeError: ...__call__() ...
        """
        if isinstance(query, str):
            if re.match('^A[0-9]{6}$', query):
                return self.find_by_id(query)
            else:
                return self.find_by_description(query, max_results, first_result)
        elif isinstance(query, (int, Integer)):
            return self.find_by_id(query)
        elif isinstance(query, (list, tuple)):
            return self.find_by_subsequence(query, max_results, first_result)

    def __repr__(self) -> str:
        r"""
        Return the representation of ``self``.

        TESTS::

            sage: oeis
            The On-Line Encyclopedia of Integer Sequences (https://oeis.org/)
        """
        return "The On-Line Encyclopedia of Integer Sequences (%s)" % oeis_url

    def find_by_id(self, ident, fetch=False):
        r"""

        INPUT:

        - ``ident`` -- a string representing the A-number of the sequence
          or an integer representing its number.

        - ``fetch`` -- (bool, default: ``False``) whether to force fetching the
          content of the sequence on the internet.

        OUTPUT:

        - The OEIS sequence whose A-number or number corresponds to ``ident``.

        EXAMPLES::

            sage: oeis.find_by_id('A000040')            # optional -- internet
            A000040: The prime numbers.

            sage: oeis.find_by_id(40)                   # optional -- internet
            A000040: The prime numbers.
        """
        sequence = OEISSequence(ident=ident)
        if fetch:
            sequence.online_update()
        return sequence

    def find_by_entry(self, entry):
        r"""

        INPUT:

        - ``entry`` -- a string corresponding to an entry in the internal format
          of the OEIS.

        OUTPUT:

        - The corresponding OEIS sequence.

        EXAMPLES::

            sage: entry = '%I A262002\n%N A262002 L.g.f.: log( Sum_{n>=0} x^n/n! * Product_{k=1..n} (k^2 + 1) ).\n%K A262002 nonn'
            sage: s = oeis.find_by_entry(entry)
            sage: s
            A262002: L.g.f.: log( Sum_{n>=0} x^n/n! * Product_{k=1..n} (k^2 + 1) ).
        """
        ident = entry[3:10]
        sequence = OEISSequence(ident=ident)
        sequence._raw = entry
        return sequence

    def find_by_description(self, description, max_results=3, first_result=0):
        r"""
        Search for OEIS sequences corresponding to the description.

        INPUT:

        - ``description`` -- (string) the description the searched sequences.

        - ``max_results`` -- (integer, default: 3) the maximum number of results
          we want. In any case, the on-line encyclopedia will not return more
          than 100 results.

        - ``first_result`` -- (integer, default: 0) allow to skip the
          ``first_result`` first results in the search, to go further.
          This is useful if you are looking for a sequence that may appear
          after the 100 first found sequences.

        OUTPUT:

        - a tuple (with fancy formatting) of at most ``max_results`` OEIS
          sequences. Those sequences can be used without the need to fetch the
          database again.

        EXAMPLES::

            sage: oeis.find_by_description('prime gap factorization') # optional -- internet
            0: A...: ...
            1: A...: ...
            2: A...: ...

            sage: prime_gaps = _[2] ; prime_gaps        # optional -- internet
            A073490: Number of prime gaps in factorization of n.

            sage: oeis('beaver')                        # optional -- internet
            0: A...: ...eaver...
            1: A...: ...eaver...
            2: A...: ...eaver...

            sage: oeis('beaver', max_results=4, first_result=2)     # optional -- internet
            0: A...: ...eaver...
            1: A...: ...eaver...
            2: A...: ...eaver...
            3: A...: ...eaver...
        """
        options = {'q': description,
                   'n': str(max_results),
                   'fmt': 'text',
                   'start': str(first_result)}
        url = oeis_url + "search?" + urlencode(options)
        sequence_list = _fetch(url).split('\n\n')[2:-1]
        return FancyTuple([self.find_by_entry(entry=_) for _ in sequence_list])

    def find_by_subsequence(self, subsequence, max_results=3, first_result=0):
        r"""
        Search for OEIS sequences containing the given subsequence.

        INPUT:

        - ``subsequence`` -- a list or tuple of integers.

        - ``max_results`` -- (integer, default: 3), the maximum of results requested.

        - ``first_result`` -- (integer, default: 0) allow to skip the
          ``first_result`` first results in the search, to go further.
          This is useful if you are looking for a sequence that may appear
          after the 100 first found sequences.

        OUTPUT:

        - a tuple (with fancy formatting) of at most ``max_results`` OEIS
          sequences. Those sequences can be used without the need to fetch the
          database again.

        EXAMPLES::

            sage: oeis.find_by_subsequence([2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377])  # optional -- internet
            0: A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
            1: A212804: Expansion of (1 - x)/(1 - x - x^2).
            2: A020695: Pisot sequence E(2,3).

            sage: fibo = _[0] ; fibo                    # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
        """
        subsequence = str(subsequence)[1:-1]
        return self.find_by_description(subsequence, max_results, first_result)

    def browse(self):
        r"""
        Open the OEIS web page in a browser.

        EXAMPLES::

            sage: oeis.browse()                         # optional -- webbrowser
        """
        import webbrowser
        webbrowser.open(oeis_url)

    def _imaginary_entry(self, ident='A999999', keywords=''):
        r"""
        This is an imaginary entry of an OEIS sequence for offline tests.

        INPUT:

        - ``ident`` -- a string representing the A-number of the sequence.

        - ``keywords`` -- a string corresponding to the keyword field of the
          sequence.

        OUTPUT:

        - a string representing the entry of the sequence.

        TESTS::

            sage: oeis._imaginary_entry().split('\n')[0]
            '%I A999999 M9999 N9999'

            sage: keywords = 'simon,cussonet'
            sage: s = oeis.find_by_entry(entry=oeis._imaginary_entry(ident='A999998', keywords=keywords))
            sage: ','.join(s.keywords()) == keywords
            True

        """
        return ('%I ' + ident + ' M9999 N9999\n'
                '%S ' + ident + ' 1,1,1,1,2,1,1,1,\n'
                '%T ' + ident + ' 1,1,1,1,1,1,1,1,1,\n'
                '%U ' + ident + ' 1,1,1,1,1,1,1,1,1\n'
                '%N ' + ident + ' The characteristic sequence of 42 plus one, starting from 38.\n'
                '%D ' + ident + ' Lewis Carroll, Alice\'s Adventures in Wonderland.\n'
                '%D ' + ident + ' Lewis Carroll, The Hunting of the Snark.\n'
                '%D ' + ident + ' Deep Thought, The Answer to the Ultimate Question of Life, The Universe, and Everything.\n'
                '%H ' + ident + ' Wikipedia, <a href="https://en.wikipedia.org/wiki/42_(number)">42 (number)</a>\n'
                '%H ' + ident + ' See. also <a href="https://trac.sagemath.org/sage_trac/ticket/42">trac ticket #42</a>\n'
                '%H ' + ident + ' Do not confuse with the sequence <a href="/A000042">A000042</a> or the sequence <a href="/A000024">A000024</a>\n'
                '%H ' + ident + ' The string http://42.com is not a link.\n'
                '%F ' + ident + ' For n big enough, s(n+1) - s(n) = 0.\n'
                '%Y ' + ident + ' Related sequences are A000042 and its friend A000024.\n'
                '%A ' + ident + ' Anonymous.\n'
                '%O ' + ident + ' 38,4\n'
                '%E ' + ident + ' This sequence does not contain errors.\n'
                '%e ' + ident + ' s(42) + s(43) = 0.\n'
                '%p ' + ident + ' Do not even try, Maple is not able to produce such a sequence.\n'
                '%t ' + ident + ' Mathematica neither.\n'
                '%o ' + ident + ' (Python)\n'
                '%o ' + ident + ' def ' + ident + '(n):\n'
                '%o ' + ident + '     assert(isinstance(n, (int, Integer))), "n must be an integer."\n'
                '%o ' + ident + '     if n < 38:\n'
                '%o ' + ident + '         raise ValueError("The value %s is not accepted." %str(n))\n'
                '%o ' + ident + '     elif n == 42:\n'
                '%o ' + ident + '         return 2\n'
                '%o ' + ident + '     else:\n'
                '%o ' + ident + '         return 1\n'
                '%K ' + ident + ' ' + keywords + '\n'
                '%C ' + ident + ' 42 is the product of the first 4 prime numbers, except 5 and perhaps 1.\n'
                '%C ' + ident + ' Apart from that, i have no comment.')

    def _imaginary_sequence(self, ident='A999999', keywords='sign,easy'):
        r"""
        This is the OEIS sequence corresponding to the imaginary entry.
        Its main purpose is to allow offline doctesting.

        INPUT:

        - ``ident`` -- a string representing the A-number of the sequence.

        - ``keywords`` -- string (default: 'sign,easy'), a list of words
          separated by commas.

        OUTPUT:

        - OEIS sequence.

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s
            A999999: The characteristic sequence of 42 plus one, starting from 38.
            sage: s[4]
            2
            sage: s(42)
            2
        """
        return self.find_by_entry(entry=self._imaginary_entry(ident=ident, keywords=keywords))


class OEISSequence(SageObject, UniqueRepresentation):
    r"""
    The class of OEIS sequences.

    This class implements OEIS sequences. They are usually produced by calls to
    the On-Line Encyclopedia of Integer Sequences, represented by the class
    :class:`OEIS`.

    .. NOTE::

        Since some sequences do not start with index 0, there is a difference
        between calling and getting item, see :meth:`__call__` for more details
        ::

            sage: sfibo = oeis('A039834')
            sage: sfibo.first_terms()[:10]              # optional -- internet
            (1, 1, 0, 1, -1, 2, -3, 5, -8, 13)

            sage: sfibo(-2)                             # optional -- internet
            1
            sage: sfibo(3)                              # optional -- internet
            2
            sage: sfibo.offsets()                       # optional -- internet
            (-2, 6)

            sage: sfibo[0]                              # optional -- internet
            1
            sage: sfibo[6]                              # optional -- internet
            -3

    .. automethod:: __call__
    """
    @staticmethod
    def __classcall__(cls, ident):
        r"""
        Canonicalize the ID of the sequence into a A-number.

        TESTS::

            sage: oeis(45) is oeis('A000045')
            True
        """
        if not isinstance(ident, str):
            ident = str(ident)
            ident = 'A000000'[:-len(ident)] + ident
        return super(OEISSequence, cls).__classcall__(cls, ident)

    def __init__(self, ident):
        r"""
        Initialize an OEIS sequence.

        There is no fetching of additional information about the sequence at
        this point, only the A-number is required to construct a sequence.

        INPUT:

        - ``ident`` -- a string representing the A-number of the sequence or an
          integer representing its number.

        TESTS::

            sage: sfibo = oeis('A039834')

            sage: s = oeis._imaginary_sequence()
        """
        self._id = ident

    def online_update(self):
        r"""
        Fetch the online OEIS to update the informations about this sequence.

        TESTS::

            sage: s = oeis._imaginary_sequence(ident='A004238')
            sage: s
            A004238: The characteristic sequence of 42 plus one, starting from 38.
            sage: s.online_update()                     # optional -- internet
            sage: s                                     # optional -- internet
            A004238: a(n) = 100*log(n) rounded to nearest integer.
        """
        options = {'q': self._id, 'n': '1', 'fmt': 'text'}
        url = oeis_url + "search?" + urlencode(options)
        self._raw = _fetch(url).split('\n\n')[2]
        try:
            del self._fields
        except AttributeError:
            pass

    def _field(self, key):
        r"""
        Return the ``key`` field of the entry of ``self``.

        This method allows to handle the ``_fields`` dictionary in a lazy way.

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s._field('C')[0]
            '42 is the product of the first 4 prime numbers, except 5 and perhaps 1.'
        """
        try:
            return self._fields[key]
        except AttributeError:
            fields = defaultdict(list)
            for line in self.raw_entry().splitlines():
                fields[line[1]].append(line[11:])
            self._fields = fields
            self.is_dead(warn_only=True)
            return self._fields[key]

    def id(self, format='A'):
        r"""
        The ID of the sequence ``self`` is the A-number that identifies
        ``self``.

        INPUT:

        - ``format`` -- (string, default: 'A').

        OUTPUT:

        - if ``format`` is set to 'A', returns a string of the form 'A000123'.
        - if ``format`` is set to 'int' returns an integer of the form 123.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.id()                                # optional -- internet
            'A000045'

            sage: f.id(format='int')                    # optional -- internet
            45

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.id()
            'A999999'
            sage: s.id(format='int')
            999999
        """
        if format == 'A':
            return self._id
        elif format == 'int':
            return int(self._id[1:].lstrip("0"))

    def __hash__(self):
        r"""
        Return the hash of ``self``, which is its numerical OEIS ID.

        This method allows unique representation of OEIS sequences.

        OUTPUT:

        - Python `int`.

        EXAMPLES::

            sage: s = oeis([1,2,3,5,8,13])[0]           # optional -- internet
            sage: hash(s)                               # optional -- internet
            45

        We have unique representation::

            sage: t = oeis(45)                          # optional -- internet
            sage: s is t                                # optional -- internet
            True
            sage: s == t                                # optional -- internet
            True

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s is oeis._imaginary_sequence()
            True
            sage: s == oeis._imaginary_sequence()
            True
        """
        return self.id(format='int')

    def raw_entry(self):
        r"""
        Return the raw entry of the sequence ``self``, in the OEIS format.

        The raw entry is fetched online if needed.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: print(f.raw_entry())                  # optional -- internet
            %I A000045 M0692 N0256
            %S A000045 0,1,1,2,3,5,8,13,21,34,55,89,144,...
            %T A000045 10946,17711,28657,46368,...
            ...

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.raw_entry() == oeis._imaginary_entry(keywords='sign,easy')
            True
        """
        try:
            return self._raw
        except AttributeError:
            self.online_update()
            return self._raw

    def name(self) -> str:
        r"""
        Return the name of the sequence ``self``.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.name()                              # optional -- internet
            'Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.name()
            'The characteristic sequence of 42 plus one, starting from 38.'
        """
        return self._field('N')[0]

    def old_IDs(self):
        r"""
        Return the IDs of the sequence ``self`` corresponding to ancestors of OEIS.

        OUTPUT:

        - a tuple of at most two strings. When the string starts with `M`, it
          corresponds to the ID of "The Encyclopedia of Integer Sequences" of
          1995. When the string starts with `N`, it corresponds to the ID of
          the "Handbook of Integer Sequences" of 1973.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.old_IDs()                           # optional -- internet
            ('M0692', 'N0256')

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.old_IDs()
            ('M9999', 'N9999')
        """
        return tuple(self._field('I')[0].split(' '))

    def offsets(self):
        r"""
        Return the offsets of the sequence ``self``.

        The first offset is the subscript of the first term in the sequence
        ``self``. When, the sequence represents the decimal expansion of a real
        number, it corresponds to the number of digits of its integer part.

        The second offset is the first term in the sequence ``self`` (starting
        from 1) whose absolute value is greater than 1. This is set to 1 if all
        the terms are 0 or +-1.

        OUTPUT:

        - tuple of two elements.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.offsets()                           # optional -- internet
            (0, 4)

            sage: f.first_terms()[:4]                   # optional -- internet
            (0, 1, 1, 2)

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.offsets()
            (38, 4)
        """
        return to_tuple(self._field('O')[0])

    def author(self):
        r"""
        Return the author of the sequence in the encyclopedia.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.author()                            # optional -- internet
            '_N. J. A. Sloane_, 1964'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.author()
            'Anonymous.'
        """
        return self._field('A')[0]

    def keywords(self):
        r"""
        Return the keywords associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings.

        EXAMPLES::

            sage: f = oeis(53) ; f                      # optional -- internet
            A000053: Local stops on New York City Broadway line (IRT #1) subway.

            sage: f.keywords()                          # optional -- internet
            ('nonn', 'fini', ...)

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.keywords()
            ('sign', 'easy')

            sage: s = oeis._imaginary_sequence(ident='A999997', keywords='nonn,hard')
            sage: s.keywords()
            ('nonn', 'hard')
        """
        return tuple(self._field('K')[0].split(','))

    def natural_object(self):
        r"""
        Return the natural object associated to the sequence ``self``.

        OUTPUT:

        - If the sequence ``self`` corresponds to the digits of a real
              number, returns the associated real number (as an element of
              RealLazyField()).

        - If the sequence ``self`` corresponds to the convergents of a
              continued fraction, returns the associated continued fraction.

        .. WARNING::

            This method forgets the fact that the returned sequence may not be
            complete.

        .. TODO::

            - ask OEIS to add a keyword telling whether the sequence comes from
              a power series, e.g. for https://oeis.org/A000182
            - discover other possible conversions.

        EXAMPLES::

            sage: g = oeis("A002852") ; g               # optional -- internet
            A002852: Continued fraction for Euler's constant (or Euler-Mascheroni constant) gamma.

            sage: x = g.natural_object() ; type(x)      # optional -- internet
            <class 'sage.rings.continued_fraction.ContinuedFraction_periodic'>

            sage: RDF(x) == RDF(euler_gamma)            # optional -- internet
            True

            sage: cfg = continued_fraction(euler_gamma)
            sage: x[:90] == cfg[:90]                    # optional -- internet
            True

        ::

            sage: ee = oeis('A001113') ; ee             # optional -- internet
            A001113: Decimal expansion of e.

            sage: x = ee.natural_object() ; x           # optional -- internet
            2.718281828459046?

            sage: x.parent()                            # optional -- internet
            Real Lazy Field

            sage: x == RR(e)                            # optional -- internet
            True

        ::

            sage: av = oeis('A087778') ; av             # optional -- internet
            A087778: Decimal expansion ... Avogadro...

            sage: av.natural_object()                   # optional -- internet
            6.022141000000000?e23

        ::

            sage: fib = oeis('A000045') ; fib           # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: x = fib.natural_object() ; x.universe()         # optional -- internet
            Non negative integer semiring

        ::

            sage: sfib = oeis('A039834') ; sfib         # optional -- internet
            A039834: a(n+2) = -a(n+1) + a(n) (signed Fibonacci numbers) with a(-2) = a(-1) = 1; or Fibonacci numbers (A000045) extended to negative indices.

            sage: x = sfib.natural_object() ; x.universe()    # optional -- internet
            Integer Ring

        TESTS::

            sage: s = oeis._imaginary_sequence(ident='A999996', keywords='nonn,cofr')
            sage: type(s.natural_object())
            <class 'sage.rings.continued_fraction.ContinuedFraction_periodic'>

            sage: s = oeis._imaginary_sequence(ident='A999995', keywords='nonn')
            sage: s.natural_object().universe()
            Non negative integer semiring

            sage: s = oeis._imaginary_sequence()
            sage: s.natural_object().universe()
            Integer Ring
        """
        if 'cofr' in self.keywords() and 'frac' not in self.keywords():
            from sage.rings.continued_fraction import continued_fraction
            return continued_fraction(self.first_terms())
        elif 'cons' in self.keywords():
            offset = self.offsets()[0]
            terms = self.first_terms() + tuple([0] * abs(offset))
            from sage.rings.real_lazy import RealLazyField
            return RealLazyField()('0' + ''.join(map(str, terms[:offset])) + '.' + ''.join(map(str, terms[offset:])))
        elif 'nonn' in self.keywords():
            from sage.structure.sequence import Sequence
            from sage.rings.semirings.non_negative_integer_semiring import NN
            return Sequence(self.first_terms(), NN)
        else:
            from sage.structure.sequence import Sequence
            from sage.rings.integer_ring import ZZ
            return Sequence(self.first_terms(), ZZ)

    def is_dead(self, warn_only=False):
        r"""
        Tell whether the sequence is dead (i.e. erroneous).

        INPUT:

        - warn_only -- (bool, default: ``False``), whether to warn when the
          sequence is dead instead of returning a boolean.

        EXAMPLES:

        A warn_only test is triggered as soon as some information on the
        sequence is queried::

            sage: s = oeis(17)
            sage: s                                 # optional -- internet
            doctest:warning
            ...
            RuntimeWarning: This sequence is dead: "A000017: Erroneous version of A032522."
            A000017: Erroneous version of A032522.

        TESTS::

            sage: s.is_dead()                       # optional -- internet
            True

            sage: t = oeis._imaginary_sequence()
            sage: t.is_dead()
            False

            sage: u = oeis._imaginary_sequence(ident='A999994', keywords='dead')
            sage: u
            doctest:warning
            ...
            RuntimeWarning: This sequence is dead: "A999994: The characteristic sequence of 42 plus one, starting from 38."
            A999994: The characteristic sequence of 42 plus one, starting from 38.

            sage: u.is_dead()
            True
        """
        if warn_only:
            if 'dead' in self.keywords():
                from warnings import warn
                warn('This sequence is dead: "{}: {}"'.format(self.id(), self.name()), RuntimeWarning)
        else:
            return 'dead' in self.keywords()

    def is_finite(self):
        r"""
        Tell whether the sequence is finite.

        Currently, OEIS only provides a keyword when the sequence is known to
        be finite. So, when this keyword is not there, we do not know whether
        it is infinite or not.

        OUTPUT:

        - ``True`` when the sequence is known to be finite.
        - ``Unknown`` otherwise.

        .. TODO::

            Ask OEIS for a keyword ensuring that a sequence is infinite.

        EXAMPLES::

            sage: s = oeis('A114288') ; s               # optional -- internet
            A114288: Lexicographically earliest solution of any 9 X 9 sudoku, read by rows.

            sage: s.is_finite()                         # optional -- internet
            True

        ::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.is_finite()                         # optional -- internet
            Unknown

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.is_finite()
            Unknown

            sage: s = oeis._imaginary_sequence(ident='A999993', keywords='nonn,finit')
            sage: s.is_finite()
            True

        """
        if 'finit' in self.keywords() or 'full' in self.keywords():
            return True
        else:
            return Unknown

    def is_full(self):
        r"""
        Tell whether the sequence ``self`` is full, that is, if all its
        elements are listed in ``self.first_terms()``.

        Currently, OEIS only provides a keyword when the sequence is known to
        be full. So, when this keyword is not there, we do not know whether
        some elements are missing or not.

        OUTPUT:

        - ``True`` when the sequence is known to be full.
        - ``Unknown`` otherwise.

        EXAMPLES::

            sage: s = oeis('A114288') ; s               # optional -- internet
            A114288: Lexicographically earliest solution of any 9 X 9 sudoku, read by rows.

            sage: s.is_full()                           # optional -- internet
            True

        ::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.is_full()                           # optional -- internet
            Unknown

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.is_full()
            Unknown

            sage: s = oeis._imaginary_sequence(ident='A999992', keywords='nonn,full,finit')
            sage: s.is_full()
            True
        """
        if 'full' in self.keywords():
            return True
        else:
            return Unknown

    @cached_method
    def first_terms(self, number=None):
        r"""

        INPUT:

        - ``number`` -- (integer or ``None``, default: ``None``) the number of
          terms returned (if less than the number of available terms). When set
          to None, returns all the known terms.

        OUTPUT:

        - tuple of integers.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.first_terms()[:10]                  # optional -- internet
            (0, 1, 1, 2, 3, 5, 8, 13, 21, 34)

        Handle dead sequences, see :trac:`17330` ::

            sage: oeis(5000).first_terms(12)              # optional -- internet
            doctest:warning
            ...
            RuntimeWarning: This sequence is dead: "A005000: Erroneous version of A006505."
            (1, 0, 0, 1, 1, 1, 11, 36, 92, 491, 2537)

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.first_terms()
            (1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
            sage: s.first_terms(5)
            (1, 1, 1, 1, 2)
        """
        fields = ['S', 'T', 'U']
        return to_tuple(" ".join(flatten([self._field(a) for a in fields])))[:number]

    def _repr_(self):
        r"""
        Print the sequence number and a short summary of this sequence.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: f = oeis(45)                          # optional -- internet
            sage: f                                     # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s
            A999999: The characteristic sequence of 42 plus one, starting from 38.
        """
        return "%s: %s" % (self.id(), self.name())

    def __call__(self, k):
        r"""
        Return the element of the sequence ``self`` with index ``k``.

        INPUT:

        - ``k`` -- integer.

        OUTPUT:

        - integer.

        .. NOTE::

            The first index of the sequence ``self`` is not necessarily zero,
            it depends on the first offset of ``self``. If the sequence
            represents the decimal expansion of a real number, the index 0
            corresponds to the digit right after the decimal point.

        EXAMPLES::

            sage: f = oeis(45)                          # optional -- internet
            sage: f.first_terms()[:10]                  # optional -- internet
            (0, 1, 1, 2, 3, 5, 8, 13, 21, 34)

            sage: f(4)                                  # optional -- internet
            3

        ::

            sage: sfibo = oeis('A039834')               # optional -- internet
            sage: sfibo.first_terms()[:10]              # optional -- internet
            (1, 1, 0, 1, -1, 2, -3, 5, -8, 13)

            sage: sfibo(-2)                             # optional -- internet
            1
            sage: sfibo(4)                              # optional -- internet
            -3
            sage: sfibo.offsets()                       # optional -- internet
            (-2, 6)

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s(38)
            1
            sage: s(42)
            2
            sage: s(2)
            Traceback (most recent call last):
            ...
            ValueError: Sequence A999999 is not defined (or known) for index 2
        """
        offset = self.offsets()[0]
        if 'cons' in self.keywords():
            offset = - offset
        n = k - offset
        if not 0 <= n < len(self.first_terms()):
            raise ValueError("Sequence %s is not defined (or known) for index %s" % (self.id(), k))
        return self.first_terms()[n]

    def __getitem__(self, i):
        r"""
        Return the ``i``th element of sequence ``self``, viewed as a tuple.

        The first element appearing in the sequence ``self``corresponds to
        ``self[0]``. Do not confuse with calling ``self(k)``.

        INPUT:

        - ``i`` -- integer.

        OUTPUT:

        - integer.

        EXAMPLES::

            sage: sfibo = oeis('A039834')               # optional -- internet
            sage: sfibo[8]                              # optional -- internet
            -8
            sage: sfibo(8)                              # optional -- internet
            -21

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s[2]
            1
            sage: s[4]
            2
            sage: s[38]
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.first_terms()[i]

    def __iter__(self):
        r"""
        Iterate over the first terms of ``self``, and raise an error if
        those first terms are exhausted and the real associated sequence
        still have terms to produce.

        OUTPUT:

        - integer.

        EXAMPLES::

            sage: p = oeis('A085823') ; p               # optional -- internet
            A085823: Numbers in which all substrings are primes.

            sage: for i in p:                           # optional -- internet
            ....:     print(i)
            2
            3
            5
            7
            23
            37
            53
            73
            373

        ::

            sage: w = oeis(7540) ; w                    # optional -- internet
            A007540: Wilson primes: primes p such that (p-1)! == -1 (mod p^2).

            sage: i = w.__iter__()                      # optional -- internet
            sage: next(i)                               # optional -- internet
            5
            sage: next(i)                               # optional -- internet
            13
            sage: next(i)                               # optional -- internet
            563
            sage: next(i)                               # optional -- internet
            Traceback (most recent call last):
            ...
            LookupError: Future values not provided by OEIS.

        ::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: for i in f:                           # optional -- internet
            ....:     print(i)
            Traceback (most recent call last):
            ...
            LookupError: Future values not provided by OEIS.

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: for i in s:
            ....:     pass
            Traceback (most recent call last):
            ...
            LookupError: Future values not provided by OEIS.

            sage: for i in s:
            ....:     if i == 2:
            ....:         print(i)
            ....:         break
            2

            sage: s = oeis._imaginary_sequence(ident='A999991', keywords='sign,full')
            sage: for i in s: pass
        """
        for x in self.first_terms():
            yield x
        if not self.is_full() is True:
            raise LookupError("Future values not provided by OEIS.")

    def references(self):
        r"""
        Return a tuple of references associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: w = oeis(7540) ; w                    # optional -- internet
            A007540: Wilson primes: primes p such that (p-1)! == -1 (mod p^2).

            sage: w.references()                        # optional -- internet
            ...Albert H. Beiler, Recreations in the Theory of Numbers, Dover, NY, 1964, p. 52.
            ...

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.references()[1]
            'Lewis Carroll, The Hunting of the Snark.'
        """
        return FancyTuple(self._field('D'))

    def links(self, browse=None, format='guess'):
        r"""
        Return, display or browse links associated to the sequence ``self``.

        INPUT:

        - ``browse`` -- an integer, a list of integers, or the word 'all'
          (default: ``None``) : which links to open in a web browser.

        - ``format`` -- string (default: 'guess') : how to display the links.

        OUTPUT:

        - tuple of strings (with fancy formatting):
            - if ``format`` is ``url``, returns a tuple of absolute links without description.
            - if ``format`` is ``html``, returns nothing but prints a tuple of clickable absolute links in their context.
            - if ``format`` is ``guess``, adapts the output to the context (command line or notebook).
            - if ``format`` is ``raw``, the links as they appear in the database, relative links are not made absolute.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.links(format='url')                             # optional -- internet
            0: https://oeis.org/A000045/b000045.txt
            1: ...
            2: ...

            sage: f.links(format='raw')                 # optional -- internet
            0: N. J. A. Sloane, <a href="/A000045/b000045.txt">The first 2000 Fibonacci numbers: Table of n, F(n) for n = 0..2000</a>
            1: ...
            2: ...

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.links(format='raw')[2]
            'Do not confuse with the sequence <a href="/A000042">A000042</a> or the sequence <a href="/A000024">A000024</a>'

            sage: s.links(format='url')[3]
            'https://oeis.org/A000024'

            sage: HTML = s.links(format="html");  HTML
            0: Wikipedia, <a href="https://en.wikipedia.org/wiki/42_(number)">42 (number)</a>
            1: See. also <a href="https://trac.sagemath.org/sage_trac/ticket/42">trac ticket #42</a>
            ...
            sage: type(HTML)
            <class 'sage.misc.html.HtmlFragment'>
        """
        def url_absolute(s):
            return re.sub(r'\"\/', '\"' + oeis_url, s)
        if browse is None:
            if format == 'guess':
                return self.links(format='url')
            elif format == 'raw':
                return FancyTuple(self._field('H'))
            elif format == 'html':
                return HtmlFragment(FancyTuple([url_absolute(_) for _ in self._field('H')]))
            elif format == 'url':
                url_list = flatten([_urls(url_absolute(string)) for string in self._field('H')])
                return FancyTuple(url_list)
        else:
            import webbrowser
            url_list = flatten([_urls(url_absolute(string)) for string in self._field('H')])
            if isinstance(browse, (int, Integer)):
                webbrowser.open(url_list[browse])
            elif isinstance(browse, (list, tuple)):
                for url_number in browse:
                    webbrowser.open(url_list[url_number])
            elif browse == 'all':
                for url in url_list:
                    webbrowser.open(url)

    def formulas(self):
        r"""
        Return a tuple of formulas associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.formulas()[2]                       # optional -- internet
            'F(n) = ((1+sqrt(5))^n - (1-sqrt(5))^n)/(2^n*sqrt(5)).'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.formulas()
            0: For n big enough, s(n+1) - s(n) = 0.

        """
        return FancyTuple(self._field('F'))

    def cross_references(self, fetch=False):
        r"""
        Return a tuple of cross references associated to the sequence
        ``self``.

        INPUT:

        - ``fetch`` -- boolean (default: ``False``).

        OUTPUT:

        - if ``fetch`` is ``False``, return a list of OEIS IDs (strings).
        - if ``fetch`` if ``True``, return a tuple of OEIS sequences.

        EXAMPLES::

            sage: nbalanced = oeis("A005598") ; nbalanced   # optional -- internet
            A005598: a(n) = 1 + Sum_{i=1..n} (n-i+1)*phi(i).

            sage: nbalanced.cross_references()              # optional -- internet
            ('A049703', 'A049695', 'A103116', 'A000010')

            sage: nbalanced.cross_references(fetch=True)    # optional -- internet
            0: A049703: a(0) = 0; for n>0, a(n) = A005598(n)/2.
            1: A049695: Array T read by diagonals; ...
            2: A103116: a(n) = Sum_{i=1..n} (n-i+1)*phi(i).
            3: A000010: Euler totient function phi(n): count numbers <= n and prime to n.

            sage: phi = _[3]                                # optional -- internet

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.cross_references()
            ('A000042', 'A000024')
        """
        ref_list = re.findall('A[0-9]{6}', " ".join(self._field('Y')))
        if fetch:
            return FancyTuple([oeis.find_by_id(_) for _ in ref_list])
        else:
            return tuple(ref_list)

    def extensions_or_errors(self):
        r"""
        Return a tuple of extensions or errors associated to the
        sequence ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: sfibo = oeis('A039834') ; sfibo       # optional -- internet
            A039834: a(n+2) = -a(n+1) + a(n) (signed Fibonacci numbers) with a(-2) = a(-1) = 1; or Fibonacci numbers (A000045) extended to negative indices.

            sage: sfibo.extensions_or_errors()[0]       # optional -- internet
            'Signs corrected by _Len Smiley_ and _N. J. A. Sloane_'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.extensions_or_errors()
            0: This sequence does not contain errors.

        """
        return FancyTuple(self._field('E'))

    def examples(self):
        r"""
        Return a tuple of examples associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: c = oeis(1203) ; c                    # optional -- internet
            A001203: Simple continued fraction expansion of Pi.

            sage: c.examples()                          # optional -- internet
            0: Pi = 3.1415926535897932384...
            1:    = 3 + 1/(7 + 1/(15 + 1/(1 + 1/(292 + ...))))
            2:    = [a_0; a_1, a_2, a_3, ...] = [3; 7, 15, 1, 292, ...].

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.examples()
            0: s(42) + s(43) = 0.
        """
        return FancyTuple(self._field('e'))

    def comments(self):
        r"""
        Return a tuple of comments associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.comments()[:8]                      # optional -- internet
            0: D. E. Knuth writes...
            1: In keeping with historical accounts...
            2: Susantha Goonatilake writes...
            3: Also sometimes called Hemachandra numbers.
            4: Also sometimes called Lamé's sequence.
            5: ...
            6: F(n+2) = number of binary sequences of length n that have no consecutive 0's.
            7: F(n+2) = number of subsets of {1,2,...,n} that contain no consecutive integers.

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.comments()
            0: 42 is the product of the first 4 prime numbers, except 5 and perhaps 1.
            1: Apart from that, i have no comment.
        """
        return FancyTuple(self._field('C'))

    def url(self):
        r"""
        Return the URL of the page associated to the sequence ``self``.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.url()                               # optional -- internet
            'https://oeis.org/A000045'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.url()
            'https://oeis.org/A999999'
        """
        return oeis_url + self.id()

    def browse(self):
        r"""
        Open the OEIS web page associated to the sequence ``self`` in a browser.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet webbrowser
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.browse()                            # optional -- internet webbrowser

        TESTS::

            sage: s = oeis._imaginary_sequence()        # optional -- webbrowser
            sage: s.browse()                            # optional -- webbrowser
        """
        import webbrowser
        webbrowser.open(self.url())

    def show(self):
        r"""
        Display most available informations about the sequence ``self``.

        EXAMPLES::

            sage: s = oeis(12345)                       # optional -- internet
            sage: s.show()                              # optional -- internet
            ID
            A012345
            <BLANKLINE>
            NAME
            Coefficients in the expansion sinh(arcsin(x)*arcsin(x)) = 2*x^2/2!+8*x^4/4!+248*x^6/6!+11328*x^8/8!+...
            <BLANKLINE>
            FIRST TERMS
            (2, 8, 248, 11328, 849312, 94857600, 14819214720, 3091936512000, 831657655349760, 280473756197529600, 115967597965430077440, 57712257892456911912960, 34039765801079493369569280)
            <BLANKLINE>
            LINKS
            0: https://oeis.org/A012345/b012345.txt
            <BLANKLINE>
            FORMULAS
            ...
            OFFSETS
            (0, 1)
            <BLANKLINE>
            URL
            https://oeis.org/A012345
            <BLANKLINE>
            AUTHOR
            Patrick Demichel (patrick.demichel(AT)hp.com)
            <BLANKLINE>

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.show()
            ID
            A999999
            <BLANKLINE>
            NAME
            The characteristic sequence of 42 plus ...
            FIRST TERMS
            (1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
            <BLANKLINE>
            COMMENTS
            0: 42 is the product of the first 4 prime numbers, except ...
            1: Apart from that, i have no comment.
            ...
        """
        for s in ['id', 'name', 'first_terms', 'comments', 'references',
                  'links', 'formulas', 'examples', 'cross_references',
                  'programs', 'keywords', 'offsets', 'url', 'old_IDs',
                  'author', 'extensions_or_errors']:
            result = getattr(self, s)()
            if result != '' and result != ('',) and result != ():
                print(re.sub('_', ' ', s).upper())
                print(str(result) + '\n')

    def programs(self, language='all', preparsing=True, keep_comments=False):
        r"""
        Return programs for the sequence ``self`` in the given ``language``.

        INPUT:

        - ``language`` -- string (default: 'all'), the chosen language.
          Possible values are 'all' for the full list, or
          any language name, for example 'sage', 'maple', 'mathematica', etc.

        Some further optional input is specific to sage code treatment:

        - ``preparsing`` -- boolean (default: ``True``) whether to preparse
          sage code
        - ``keep_comments`` -- boolean (default: ``False``) whether to keep
          comments in sage code

        OUTPUT:

        If ``language`` is ``'all'``, this returns a sorted list of pairs
        (language, code), where every language can appear several times.

        Otherwise, this returns  a list of programs in the ``language``,
        each program being a tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: ee = oeis('A001113') ; ee             # optional -- internet
            A001113: Decimal expansion of e.

            sage: ee.programs('pari')[0]                # optional -- internet
            0: default(realprecision, 50080); x=exp(1); for (n=1, 50000, d=floor(x); x=(x-d)*10; write("b001113.txt", n, " ", d)); \\ _Harry J. Smith_, Apr 15 2009

            sage: G = oeis.find_by_id('A27642')   # optional -- internet
            sage: G.programs('all')               # optional -- internet
            [('haskell', ...),
             ('magma', ...),
             ...
             ('python', ...),
             ('sage', ...)]

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.programs()
            [('maple', ...),
            ('mathematica', ...),
            ('python',
            0: def A999999(n):
            1:     assert(isinstance(n, (int, Integer))), "n must be an integer."
            2:     if n < 38:
            3:         raise ValueError("The value %s is not accepted." %str(n))
            4:     elif n == 42:
            5:         return 2
            6:     else:
            7:         return 1)]

            sage: s.programs('maple')[0]
            0: Do not even try, Maple is not able to produce such a sequence.

            sage: s.programs('mathematica')[0]
            0: Mathematica neither.
        """
        language = language.lower()
        if language == "maple":
            return [FancyTuple(self._field('p'))]
        elif language == "mathematica":
            return [FancyTuple(self._field('t'))]
        if language == 'sagemath':
            language = 'sage'
        if language == 'all':
            table = (('maple', FancyTuple(self._field('p'))),
                     ('mathematica', FancyTuple(self._field('t'))))
            table = [(lang, code) for lang, code in table if code]
        else:
            table = []

        known_langs = ['sage', 'python', 'scheme']

        def is_starting_line(line):
            """
            Help to split the big OEIS code block into blocks by language.

            This returns ``None`` if ``line`` is not a starting line.
            """
            if not line.startswith('('):
                return None
            if ')' not in line:
                return None
            end = line.index(')')
            language = line[1:end].lower()  # normalise the language names
            if '(' in language:
                return None
            for special in known_langs:
                if special in language:
                    language = special
            if ' ' in language:  # get rid of language versions
                language = language.split(' ')[0]
            if language == 'c#' or language == 'c++':
                language = 'c'
            if language.isalnum():
                return (language, end)
            return None

        def filter_sage(lines):
            """
            Remove comments and preparse if required, only for sage code.

            This is an iterator.
            """
            for line in lines:
                if keep_comments or not line.strip().startswith('#'):
                    if preparsing:
                        yield preparse(line)
                    else:
                        yield line

        def flush_to_table(language, code_lines):
            """
            Put a list of code lines into the appropriate box of the table.

            With special treatment for sage code blocks.
            """
            if language == 'sage':
                table.append((language, FancyTuple(filter_sage(code_lines))))
            elif language is not None:
                table.append((language, FancyTuple(code_lines)))

        programs = FancyTuple(self._field('o'))
        code_lines = []
        old_language = None
        for line in programs:
            new_language = is_starting_line(line)
            if new_language is not None:
                # flush the stock of code lines if any
                flush_to_table(old_language, code_lines)
                # start new stock of code lines
                old_language, end = new_language
                rest = line[end + 1:].strip()
                code_lines = [rest] if rest else []
            else:
                code_lines.append(line)
        flush_to_table(old_language, code_lines)

        if language == 'all':
            return sorted(table)
        return sorted(prog for la, prog in table if la == language)

    def test_compile_sage_code(self):
        """
        Try to compile the extracted sage code, if there is any.

        If there are several sage code fields, they are all considered.

        Dead sequences are considered to compile correctly by default.

        This returns ``True`` if the code compiles, and raises an error
        otherwise.

        EXAMPLES:

        One correct sequence::

            sage: s = oeis.find_by_id('A027642')     # optional -- internet
            sage: s.test_compile_sage_code()    # optional -- internet
            True

        One dead sequence::

            sage: s = oeis.find_by_id('A000154')     # optional -- internet
            sage: s.test_compile_sage_code()    # optional -- internet
            doctest:warning
            ...
            RuntimeWarning: This sequence is dead: ...
            True
        """
        if self.is_dead():
            return True
        filt = self.programs(language='sage')
        if filt:
            for v in filt:
                tp = tmp_filename(ext='.sage')
                _ = compile('\n'.join(v), tp, 'exec')
        return True


class FancyTuple(tuple):
    r"""
    This class inherits from ``tuple``, it allows to nicely print tuples whose
    elements have a one line representation.

    EXAMPLES::

        sage: from sage.databases.oeis import FancyTuple
        sage: t = FancyTuple(['zero', 'one', 'two', 'three', 4]) ; t
        0: zero
        1: one
        2: two
        3: three
        4: 4

        sage: t[2]
        'two'
    """
    def __repr__(self):
        r"""
        Print the tuple with one value per line, where each line
        begins with the index of the value in ``self``.

        EXAMPLES::

            sage: from sage.databases.oeis import FancyTuple
            sage: t = FancyTuple(['zero', 'one', 'two', 'three', 4]) ; t
            0: zero
            1: one
            2: two
            3: three
            4: 4

            sage: t = FancyTuple(['Français', 'Español', '中文']) ; t
            0: Français
            1: Español
            2: 中文
        """
        length = len(str(len(self) - 1))
        return '\n'.join('{0:>{1}}: {2}'.format(i, length, item) for i, item in enumerate(self))

    def __getslice__(self, i, j):
        r"""
        The slice of a FancyTuple remains a FancyTuple.

        EXAMPLES::

            sage: from sage.databases.oeis import FancyTuple
            sage: t = FancyTuple(['zero', 'one', 'two', 'three', 4])
            sage: t[-2:]
            0: three
            1: 4

        TESTS::

            sage: t = ('é', 'è', 'à', 'ç')
            sage: FancyTuple(t)[2:4]
            0: à
            1: ç
        """
        return self.__getitem__(slice(i, j))

    def __getitem__(self, x):
        r"""
        If ``x`` is a slice return the corresponding sub FancyTuple,
        else return the `̀`x``-th item of ``self``.

        TESTS::

            sage: from sage.databases.oeis import FancyTuple
            sage: t = ('é', 'è', 'à', 'ç')
            sage: ft = FancyTuple(t)
            sage: ft[0] == 'é'
            True
            sage: ft[-1] == 'ç'
            True

        Check that :trac:`26997` is fixed::

            sage: FancyTuple([[1,2,3],(4,5,6)])
            0: [1, 2, 3]
            1: (4, 5, 6)
        """
        res = tuple.__getitem__(self, x)
        if isinstance(x, slice):
            res = FancyTuple(res)
        return res


oeis = OEIS()
