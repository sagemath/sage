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

EXAMPLES::

        sage: oeis
        The On-Line Encyclopedia of Integer Sequences (http://oeis.org/)

What about a sequence starting with `3, 7, 15, 1` ?

::

    sage: search = oeis([3, 7, 15, 1], max_results=4) ; search  # optional -- internet
    0: A001203: Continued fraction expansion of Pi.
    1: A082495: a(n) = (2^n - 1) mod n.
    2: A165416: Irregular array read by rows: The n-th row contains those distinct positive integers that each, when written in binary, occurs as a substring in binary n.
    3: A246674: Run Length Transform of A000225.

    sage: [u.id() for u in search]                      # optional -- internet
    ['A001203', 'A082495', 'A165416', 'A246674']
    sage: c = search[0] ; c                             # optional -- internet
    A001203: Continued fraction expansion of Pi.

::

    sage: c.first_terms(15)                             # optional -- internet
    (3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1)

    sage: c.examples()                                  # optional -- internet
    0: Pi = 3.1415926535897932384...
    1:    = 3 + 1/(7 + 1/(15 + 1/(1 + 1/(292 + ...))))
    2:    = [a_0; a_1, a_2, a_3, ...] = [3; 7, 15, 1, 292, ...]

    sage: c.comments()                                  # optional -- internet
    0: The first 5,821,569,425 terms were computed by _Eric W. Weisstein_ on Sep 18 2011.
    1: The first 10,672,905,501 terms were computed by _Eric W. Weisstein_ on Jul 17 2013.
    2: The first 15,000,000,000 terms were computed by _Eric W. Weisstein_ on Jul 27 2013.

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

    sage: p.cross_references(fetch=True)                # optional -- internet
    0: A000798: Number of different quasi-orders (or topologies, or transitive digraphs) with n labeled elements.
    1: A001035: Number of partially ordered sets ("posets") with n labeled elements (or labeled acyclic transitive digraphs).
    2: A001930: Number of topologies, or transitive digraphs with n unlabeled nodes.
    3: A006057: Number of topologies on n labeled points satisfying axioms T_0-T_4.
    4: A079263: Number of constrained mixed models with n factors.
    5: A079265: Number of antisymmetric transitive binary relations on n unlabeled points.


What does the Taylor expansion of the `e^(e^x-1)`` function have to do with
primes ?

::

    sage: x = var('x') ; f(x) = e^(e^x - 1)
    sage: L = [a*factorial(b) for a,b in taylor(f(x), x, 0, 20).coefficients()] ; L
    [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597,
    27644437, 190899322, 1382958545, 10480142147, 82864869804, 682076806159,
    5832742205057, 51724158235372]

    sage: oeis(L)                                       # optional -- internet
    0: A000110: Bell or exponential numbers: number of ways to partition a set of n labeled elements.

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

#*****************************************************************************
#       Copyright (C) 2012 Thierry Monteil <sage!lma.metelu.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function


from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.misc.misc import verbose
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.unknown import Unknown
from sage.misc.misc import embedded
from sage.misc.html import HtmlFragment
from collections import defaultdict
from urllib import urlopen, urlencode
import re

oeis_url = 'http://oeis.org/'


def _fetch(url):
    r"""
    Fetch the given ``url``.

    INPUT:

    - ``url`` - a string corresponding to the URL to be fetched.

    OUTPUT:

    - a string representing the fetched web page.

    TESTS::

        sage: from sage.databases.oeis import _fetch, oeis_url
        sage: _fetch(oeis_url + 'hints.html')[-8:-1]            # optional -- internet
        '</html>'
    """
    try:
        verbose("Fetching URL %s ..." % url, caller_name='OEIS')
        f = urlopen(url)
        result = f.read()
        f.close()
        return result
    except IOError as msg:
        raise IOError("%s\nError fetching %s." % (msg, url))


def _urls(html_string):
    r"""
    Return the list of URLs contained in ``html_string``.

    Only URLs provided by HTML hyperlinks (``href`` attribute of ``<a>`` tags)
    in are returned, not text strings starting with ``http://``.

    INPUT:

    - ``html_string`` - a string representing some HTML code.

    OUTPUT:

    - a list of (string) URLs contained in ``html_string``.

    EXAMPLES::

        sage: from sage.databases.oeis import _urls
        sage: html = 'http://example.com is not a link, but <a href="http://sagemath.org/">sagemath</a> is'
        sage: _urls(html)
        ['http://sagemath.org/']

    """
    urls = []
    from HTMLParser import HTMLParser

    class MyHTMLParser(HTMLParser):
        def handle_starttag(self, tag, attrs):
            if tag == 'a':
                for attr in attrs:
                    if attr[0] == 'href':
                        urls.append(attr[1])
    MyHTMLParser().feed(html_string)
    return urls


to_tuple = lambda string: tuple(Integer(x) for x in string.split(",") if x)


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

    - if ``query`` is a list of integers, returns a tuple of OEIS sequences
      containing it as a subsequence. Those sequences can be used without
      the need to fetch the database again.

    EXAMPLES::

        sage: oeis
        The On-Line Encyclopedia of Integer Sequences (http://oeis.org/)

    A particular sequence can be called by its A-number or number::

        sage: oeis('A000040')                           # optional -- internet
        A000040: The prime numbers.

        sage: oeis(45)                                  # optional -- internet
        A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

    The database can be searched by subsequence::

        sage: search = oeis([1,2,3,5,8,13]) ; search    # optional -- internet
        0: A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
        1: A027926: Triangular array T read by rows: T(n,0) = T(n,2n) = 1 for n >= 0; T(n,1) = 1 for n >= 1; T(n,k) = T(n-1,k-2) + T(n-1,k-1) for k = 2..2n-1, n >= 2.
        2: A001129: Iccanobif numbers: reverse digits of two previous terms and add.

        sage: fibo = search[0]                         # optional -- internet

        sage: fibo.name()                               # optional -- internet
        'Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.'

        sage: fibo.first_terms()                        # optional -- internet
        (0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
        1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393,
        196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887,
        9227465, 14930352, 24157817, 39088169)

        sage: fibo.cross_references()[0]                # optional -- internet
        'A039834'

        sage: fibo == oeis(45)                          # optional -- internet
        True

        sage: sfibo = oeis('A039834')                   # optional -- internet
        sage: sfibo.first_terms()                       # optional -- internet
        (1, 1, 0, 1, -1, 2, -3, 5, -8, 13, -21, 34, -55, 89, -144, 233,
        -377, 610, -987, 1597, -2584, 4181, -6765, 10946, -17711, 28657,
        -46368, 75025, -121393, 196418, -317811, 514229, -832040, 1346269,
        -2178309, 3524578, -5702887, 9227465, -14930352, 24157817)

        sage: sfibo.first_terms(absolute_value=True)[2:20] == fibo.first_terms()[:18]   # optional -- internet
        True

        sage: fibo.formulas()[4]                        # optional -- internet
        'F(n) = F(n-1) + F(n-2) = -(-1)^n F(-n).'

        sage: fibo.comments()[1]                        # optional -- internet
        "F(n+2) = number of binary sequences of length n that have no
        consecutive 0's."

        sage: fibo.links()[0]                           # optional -- internet
        'http://oeis.org/A000045/b000045.txt'

    The database can be searched by description::

        sage: oeis('prime gap factorization', max_results=4)                # optional -- internet
        0: A073491: Numbers having no prime gaps in their factorization.
        1: A073490: Number of prime gaps in factorization of n.
        2: A073492: Numbers having at least one prime gap in their factorization.
        3: A073493: Numbers having exactly one prime gap in their factorization.

    .. WARNING::

        The following will fetch the OEIS database twice (once for searching the
        database, and once again for creating the sequence ``fibo``)::

            sage: oeis([1,2,3,5,8,13])                  # optional -- internet
            0: A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
            1: A027926: Triangular array T read by rows: T(n,0) = T(n,2n) = 1 for n >= 0; T(n,1) = 1 for n >= 1; T(n,k) = T(n-1,k-2) + T(n-1,k-1) for k = 2..2n-1, n >= 2.
            2: A001129: Iccanobif numbers: reverse digits of two previous terms and add.

            sage: fibo = oeis('A000045')                # optional -- internet

        Do not do this, it is slow, it costs bandwidth and server resources !
        Instead, do the following, to reuse the result of the search to create
        the sequence::

            sage: oeis([1,2,3,5,8,13])                  # optional -- internet
            0: A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
            1: A027926: Triangular array T read by rows: T(n,0) = T(n,2n) = 1 for n >= 0; T(n,1) = 1 for n >= 1; T(n,k) = T(n-1,k-2) + T(n-1,k-1) for k = 2..2n-1, n >= 2.
            2: A001129: Iccanobif numbers: reverse digits of two previous terms and add.

            sage: fibo = _[0]                           # optional -- internet
    """

    def __call__(self, query, max_results=3, first_result=0):
        r"""
        See the documentation of :class:`OEIS`.

        TESTS::

            sage: oeis()
            Traceback (most recent call last):
            ...
            TypeError: __call__() takes at least 2 arguments (1 given)
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

    def __repr__(self):
        r"""
        Return the representation of ``self``.

        TESTS::

            sage: oeis
            The On-Line Encyclopedia of Integer Sequences (http://oeis.org/)
        """
        return "The On-Line Encyclopedia of Integer Sequences (%s)" % oeis_url

    def find_by_id(self, ident):
        r"""

        INPUT:

        - ``ident`` - a string representing the A-number of the sequence
          or an integer representing its number.

        OUTPUT:

        - The OEIS sequence whose A-number or number corresponds to
          ``ident``.

        EXAMPLES::

            sage: oeis.find_by_id('A000040')            # optional -- internet
            A000040: The prime numbers.

            sage: oeis.find_by_id(40)                   # optional -- internet
            A000040: The prime numbers.
        """
        if not isinstance(ident, str):
            ident = str(ident)
            ident = 'A000000'[:-len(ident)] + ident
        options = {'q': ident, 'n': '1', 'fmt': 'text'}
        url = oeis_url + "search?" + urlencode(options)
        sequence = _fetch(url).split('\n\n')[2]
        return OEISSequence(sequence)

    def find_by_description(self, description, max_results=3, first_result=0):
        r"""
        Search for OEIS sequences corresponding to the description.

        INPUT:

        - ``description`` - (string) the description the searched sequences.

        - ``max_results`` - (integer, default: 3) the maximum number of results
          we want. In any case, the on-line encyclopedia will not return more
          than 100 results.

        - ``first_result`` - (integer, default: 0) allow to skip the
          ``first_result`` first results in the search, to go further.
          This is useful if you are looking for a sequence that may appear
          after the 100 first found sequences.

        OUTPUT:

        - a tuple (with fancy formatting) of at most ``max_results`` OEIS
          sequences. Those sequences can be used without the need to fetch the
          database again.

        EXAMPLES::

            sage: oeis.find_by_description('prime gap factorization')       # optional -- internet
            0: A073491: Numbers having no prime gaps in their factorization.
            1: A073490: Number of prime gaps in factorization of n.
            2: A073492: Numbers having at least one prime gap in their factorization.

            sage: prime_gaps = _[1] ; prime_gaps        # optional -- internet
            A073490: Number of prime gaps in factorization of n.

        ::

            sage: oeis('beaver')                        # optional -- internet
            0: A028444: Busy Beaver sequence, or Rado's sigma function: ...
            1: A060843: Busy Beaver problem: a(n) = maximal number of steps ...
            2: A131956: Busy Beaver variation: maximum number of steps for ...

            sage: oeis('beaver', max_results=4, first_result=2)     # optional -- internet
            0: A131956: Busy Beaver variation: maximum number of steps for ...
            1: A141475: Number of Turing machines with n states following ...
            2: A131957: Busy Beaver sigma variation: maximum number of 1's ...
            3: A052200: Number of n-state, 2-symbol, d+ in {LEFT, RIGHT}, ...
        """
        options = {'q': description,
                   'n': str(max_results),
                   'fmt': 'text',
                   'start': str(first_result)}
        url = oeis_url + "search?" + urlencode(options)
        sequence_list = _fetch(url).split('\n\n')[2:-1]
        return FancyTuple([OEISSequence(_) for _ in sequence_list])

    def find_by_subsequence(self, subsequence, max_results=3, first_result=0):
        r"""
        Search for OEIS sequences containing the given subsequence.

        INPUT:

        - ``subsequence`` - a list of integers.

        - ``max_results`` - (integer, default: 3), the maximum of results requested.

        - ``first_result`` - (integer, default: 0) allow to skip the
          ``first_result`` first results in the search, to go further.
          This is useful if you are looking for a sequence that may appear
          after the 100 first found sequences.

        OUTPUT:

        - a tuple (with fancy formatting) of at most ``max_results`` OEIS
          sequences. Those sequences can be used without the need to fetch the
          database again.

        EXAMPLES::

            sage: oeis.find_by_subsequence([2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377]) # optional -- internet
            0: A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.
            1: A177194: Fibonacci numbers whose decimal expression does not contain any digit 0.
            2: A212804: Expansion of (1-x)/(1-x-x^2).

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

    def _imaginary_entry(self, keywords=''):
        r"""
        This is an imaginary entry of an OEIS sequence for offline tests.

        INPUT:

        - ``keywords`` - a string corresponding to the keyword field of the
          sequence.

        OUTPUT:

        - a string representing the entry of the sequence.

        TESTS::

            sage: oeis._imaginary_entry().split('\n')[0]
            '%I A999999 M9999 N9999'

            sage: from sage.databases.oeis import OEISSequence
            sage: keywords = 'simon,cussonet'
            sage: s = OEISSequence(oeis._imaginary_entry(keywords))
            sage: ','.join(s.keywords()) == keywords
            True

        """
        return ('%I A999999 M9999 N9999\n'
                '%S A999999 1,1,1,1,1,1,1,1,\n'
                '%T A999999 1,1,1,1,1,1,1,1,1,\n'
                '%U A999999 1,1,1,1,1,1,1,1,1\n'
                '%V A999999 1,1,1,1,-1,1,1,1,\n'
                '%W A999999 1,1,1,1,1,1,1,1,1,\n'
                '%X A999999 1,1,1,1,1,1,1,1,1\n'
                '%N A999999 The opposite of twice the characteristic sequence of 42 plus one, starting from 38.\n'
                '%D A999999 Lewis Carroll, Alice\'s Adventures in Wonderland.\n'
                '%D A999999 Lewis Carroll, The Hunting of the Snark.\n'
                '%D A999999 Deep Thought, The Answer to the Ultimate Question of Life, The Universe, and Everything.\n'
                '%H A999999 Wikipedia, <a href="http://en.wikipedia.org/wiki/42_(number)">42 (number)</a>\n'
                '%H A999999 See. also <a href="https://trac.sagemath.org/sage_trac/ticket/42">trac ticket #42</a>\n'
                '%H A999999 Do not confuse with the sequence <a href="/A000042">A000042</a> or the sequence <a href="/A000024">A000024</a>\n'
                '%H A999999 The string http://42.com is not a link.\n'
                '%F A999999 For n big enough, s(n+1) - s(n) = 0.\n'
                '%Y A999999 Related sequences are A000042 and its friend A000024.\n'
                '%A A999999 Anonymous.\n'
                '%O A999999 38,4\n'
                '%E A999999 This sequence does not contain errors.\n'
                '%e A999999 s(42) + s(43) = 0.\n'
                '%p A999999 Do not even try, Maple is not able to produce such a sequence.\n'
                '%t A999999 Mathematica neither.\n'
                '%o A999999 (Python)\n'
                '%o A999999 def A999999(n):\n'
                '%o A999999     assert(isinstance(n, (int, Integer))), "n must be an integer."\n'
                '%o A999999     if n < 38:\n'
                '%o A999999         raise ValueError("The value %s is not accepted." %str(n)))\n'
                '%o A999999     elif n == 42:\n'
                '%o A999999         return -1\n'
                '%o A999999     else:\n'
                '%o A999999         return 1\n'
                '%K A999999 ' + keywords + '\n'
                '%C A999999 42 is the product of the first 4 prime numbers, except 5 and perhaps 1.\n'
                '%C A999999 Apart from that, i have no comment.')

    def _imaginary_sequence(self, keywords='sign,easy'):
        r"""
        This is the OEIS sequence corresponding to the imaginary entry.
        Its main purpose is to allow offline doctesting.

        INPUT:

        - ``keywords`` - string (default: 'sign,easy'), a list of words
          separated by commas.

        OUTPUT:

        - OEIS sequence.

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s
            A999999: The opposite of twice the characteristic sequence of 42 plus one, starting from 38.
            sage: s[4]
            -1
            sage: s(42)
            -1
        """
        return OEISSequence(self._imaginary_entry(keywords))


class OEISSequence(SageObject):
    r"""
    The class of OEIS sequences.

    This class implements OEIS sequences. Such sequences are produced from a
    string in the OEIS format. They are usually produced by calls to the
    On-Line Encyclopedia of Integer Sequences, represented by the class
    :class:`OEIS`.

    .. NOTE::

        Since some sequences do not start with index 0, there is a difference
        between calling and getting item, see :meth:`__call__` for more details
        ::

            sage: sfibo = oeis('A039834')               # optional -- internet
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

    def __init__(self, entry):
        r"""
        Initializes an OEIS sequence.

        TESTS::

            sage: sfibo = oeis('A039834')               # optional -- internet

        Handle dead sequences: see  :trac:`17330` ::

            sage: oeis(17)                              # optional -- internet
            .. RuntimeWarning: This sequence is dead  "A000017: Erroneous version of A032522."
            A000017: Erroneous version of A032522.

            sage: s = oeis._imaginary_sequence()
        """
        self._raw = entry
        self._id = entry[3:10]
        self._fields = defaultdict(list)
        for line in entry.splitlines():
            self._fields[line[1]].append(line[11:])
        if 'dead' in self.keywords(): 
            ("This sequence is dead: \""+self.name()+"\"")
            from warnings import warn
            warn('This sequence is dead  "'+self.id()+": "+self.name()+'"', RuntimeWarning)


    def id(self, format='A'):
        r"""
        The ID of the sequence ``self`` is the A-number that identifies
        ``self``.

        INPUT:

        - ``format`` - (string, default: 'A').

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
            return Integer(self._id[1:].lstrip("0"))

    def raw_entry(self):
        r"""
        Return the raw entry of the sequence ``self``, in the OEIS format.

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
            sage: s.raw_entry() == oeis._imaginary_entry('sign,easy')
            True
        """
        return self._raw

    def name(self):
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
            'The opposite of twice the characteristic sequence of 42 plus one, starting from 38.'
        """
        return self._fields['N'][0]

    def old_IDs(self):
        r"""
        Returns the IDs of the sequence ``self`` corresponding to ancestors of OEIS.

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
        return tuple(self._fields['I'][0].split(' '))

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
        return to_tuple(self._fields['O'][0])

    def author(self):
        r"""
        Returns the author of the sequence in the encyclopedia.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.author()                            # optional -- internet
            '_N. J. A. Sloane_, Apr 30 1991'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.author()
            'Anonymous.'
        """
        return self._fields['A'][0]

    def keywords(self):
        r"""
        Return the keywords associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings.

        EXAMPLES::

            sage: f = oeis(53) ; f                      # optional -- internet
            A000053: Local stops on New York City Broadway line (IRT #1) subway.

            sage: f.keywords()                          # optional -- internet
            ('nonn', 'fini', 'full')

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.keywords()
            ('sign', 'easy')

            sage: s = oeis._imaginary_sequence(keywords='nonn,hard')
            sage: s.keywords()
            ('nonn', 'hard')
        """
        return tuple(self._fields['K'][0].split(','))

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
              a power series, e.g. for http://oeis.org/A000182
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
            A087778: Decimal expansion of Avogadro's constant.

            sage: av.natural_object()                   # optional -- internet
            6.022141000000000?e23

        ::

            sage: fib = oeis('A000045') ; fib           # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: x = fib.natural_object() ; x.universe()         # optional -- internet
            Non negative integer semiring

        ::

            sage: sfib = oeis('A039834') ; sfib         # optional -- internet
            A039834: a(n+2) = -a(n+1)+a(n) (signed Fibonacci numbers); or Fibonacci numbers (A000045) extended to negative indices.

            sage: x = sfib.natural_object() ; x.universe()    # optional -- internet
            Integer Ring

        TESTS::

            sage: s = oeis._imaginary_sequence('nonn,cofr')
            sage: type(s.natural_object())
            <class 'sage.rings.continued_fraction.ContinuedFraction_periodic'>

            sage: s = oeis._imaginary_sequence('nonn')
            sage: s.natural_object().universe()
            Non negative integer semiring

            sage: s = oeis._imaginary_sequence()
            sage: s.natural_object().universe()
            Integer Ring
        """
        if 'cofr' in self.keywords() and not 'frac' in self.keywords():
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

    def is_finite(self):
        r"""
        Tells whether the sequence is finite.

        Currently, OEIS only provides a keyword when the sequence is known to
        be finite. So, when this keyword is not there, we do not know whether
        it is infinite or not.

        OUTPUT:

        - Returns ``True`` when the sequence is known to be finite.
        - Returns ``Unknown`` otherwise.

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

            sage: s = oeis._imaginary_sequence('nonn,finit')
            sage: s.is_finite()
            True

        """
        if 'finit' in self.keywords() or 'full' in self.keywords():
            return True
        else:
            return Unknown

    def is_full(self):
        r"""
        Tells whether the sequence ``self`` is full, that is, if all its
        elements are listed in ``self.first_terms()``.

        Currently, OEIS only provides a keyword when the sequence is known to
        be full. So, when this keyword is not there, we do not know whether
        some elements are missing or not.

        OUTPUT:

        - Returns ``True`` when the sequence is known to be full.
        - Returns ``Unknown`` otherwise.

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

            sage: s = oeis._imaginary_sequence('nonn,full,finit')
            sage: s.is_full()
            True
        """
        if 'full' in self.keywords():
            return True
        else:
            return Unknown

    @cached_method
    def first_terms(self, number=None, absolute_value=False):
        r"""

        INPUT:

        - ``number`` - (integer or ``None``, default: ``None``) the number of
          terms returned (if less than the number of available terms). When set
          to None, returns all the known terms.

        - ``absolute_value`` - (bool, default: ``False``) when a sequence has
          negative entries, OEIS also stores the absolute values of its first
          terms, when ``absolute_value`` is set to ``True``, you will get them.

        OUTPUT:

        - tuple of integers.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.first_terms()[:10]                  # optional -- internet
            (0, 1, 1, 2, 3, 5, 8, 13, 21, 34)

        Handle dead sequences: see  :trac:`17330` ::

            sage: oeis(17).first_terms(12)              # optional -- internet  
            oeis(17).first_terms(12)
            .. RuntimeWarning: This sequence is dead  "A000017: Erroneous version of A032522."
            warn('This sequence is dead  "'+self.id()+": "+self.name()+'"', RuntimeWarning)
            (1, 0, 0, 2, 2, 4, 8, 4, 16, 12, 48, 80)            

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.first_terms()
            (1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
            sage: s.first_terms(5)
            (1, 1, 1, 1, -1)
            sage: s.first_terms(5, absolute_value=True)
            (1, 1, 1, 1, 1)

            sage: s = oeis._imaginary_sequence(keywords='full')
            sage: s(40)
            Traceback (most recent call last):
            ...
            TypeError: You found a sign inconsistency, please contact OEIS

            sage: s = oeis._imaginary_sequence(keywords='sign,full')
            sage: s(40)
            1

            sage: s = oeis._imaginary_sequence(keywords='nonn,full')
            sage: s(42)
            1
        """
        if absolute_value or ('nonn' in self.keywords()) or ('dead' in self.keywords()):
            fields = ['S', 'T', 'U']
        elif ('sign' in self.keywords()):
            fields = ['V', 'W', 'X']
        else:
            raise TypeError("You found a sign inconsistency, please contact OEIS")
        return to_tuple(" ".join(flatten([self._fields[a] for a in fields])))[:number]

    def _repr_(self):
        r"""
        Prints the sequence number and a short summary of this sequence.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: f = oeis(45)                          # optional -- internet
            sage: f                                     # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s
            A999999: The opposite of twice the characteristic sequence of 42 plus one, starting from 38.
        """
        return "%s: %s" % (self.id(), self.name())

    def __call__(self, k):
        r"""
        Return the element of the sequence ``self`` with index ``k``.

        INPUT:

        - ``k`` - integer.

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
            -1
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

        - ``i`` - integer.

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
            -1
            sage: s[38]
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.first_terms()[i]

    def __iter__(self):
        r"""
        Iterates over the first terms of ``self``, and raises an error if
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
            ....:     if i == -1:
            ....:         print(i)
            ....:         break
            -1

            sage: s = oeis._imaginary_sequence(keywords='sign,full')
            sage: for i in s: pass
        """
        for x in self.first_terms():
            yield x
        if not self.is_full():
            raise LookupError("Future values not provided by OEIS.")

    def __eq__(self, other):
        r"""
        Returns ``True`` if ``self`` is equal to ``other`` and ``False``
        otherwise.  Two integer sequences are considered equal if they have the
        same OEIS ID.

        INPUT:

        - ``other`` - an oeis sequence.

        OUTPUT:

        - boolean.

        EXAMPLES::

            sage: oeis([1,2,3,5,8,13])[0] == oeis(45)   # optional -- internet
            True

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s == oeis._imaginary_sequence()
            True

        """
        return self.id() == other.id()

    def __ne__(self, other):
        r"""
        Returns ``True`` if ``self`` has a different OEIS ID than ``other`` and
        ``False`` otherwise.

        INPUT:

        - ``other`` - an oeis sequence.

        OUTPUT:

        - boolean.

        EXAMPLES::

            sage: oeis([1,2,3,5,8,13])[0] != oeis(40)   # optional -- internet
            True

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s != oeis._imaginary_sequence()
            False
        """
        return not self == other

    def references(self):
        r"""
        Return a tuple of references associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: w = oeis(7540) ; w                    # optional -- internet
            A007540: Wilson primes: primes p such that (p-1)! == -1 (mod p^2).

            sage: w.references()                        # optional -- internet
            0: A. H. Beiler, Recreations in the Theory of Numbers, Dover, NY, 1964, p. 52.
            1: C. Clawson, Mathematical Mysteries, Plenum Press, 1996, p. 180.
            2: R. Crandall and C. Pomerance, Prime Numbers: A Computational Perspective, Springer, NY, 2001; see p. 29.
            3: G. H. Hardy and E. M. Wright, An Introduction to the Theory of Numbers, 5th ed., Oxford Univ. Press, 1979, th. 80.
            ...

            sage: _[0]                                  # optional -- internet
            'A. H. Beiler, Recreations in the Theory of Numbers, Dover, NY, 1964, p. 52.'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.references()[1]
            'Lewis Carroll, The Hunting of the Snark.'
        """
        return FancyTuple(self._fields['D'])

    def links(self, browse=None, format='guess'):
        r"""
        Return, display or browse links associated to the sequence ``self``.

        INPUT:

        - ``browse`` - an integer, a list of integers, or the word 'all'
          (default: ``None``) : which links to open in a web browser.

        - ``format`` - string (default: 'guess') : how to display the links.

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
            0: http://oeis.org/A000045/b000045.txt
            1: http://www.schoolnet.ca/vp-pv/amof/e_fiboI.htm
            ...

            sage: f.links(format='raw')                 # optional -- internet
            0: N. J. A. Sloane, <a href="/A000045/b000045.txt">The first 2000 Fibonacci numbers: Table of n, F(n) for n = 0..2000</a>
            1: Amazing Mathematical Object Factory, <a href="http://www.schoolnet.ca/vp-pv/amof/e_fiboI.htm">Information on the Fibonacci sequences</a>
            ...

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.links(format='raw')[2]
            'Do not confuse with the sequence <a href="/A000042">A000042</a> or the sequence <a href="/A000024">A000024</a>'

            sage: s.links(format='url')[3]
            'http://oeis.org/A000024'

            sage: HTML = s.links(format="html");  HTML
            0: Wikipedia, <a href="http://en.wikipedia.org/wiki/42_(number)">42 (number)</a>
            1: See. also <a href="https://trac.sagemath.org/sage_trac/ticket/42">trac ticket #42</a>
            ...
            sage: type(HTML)
            <class 'sage.misc.html.HtmlFragment'>
        """
        url_absolute = lambda s: re.sub('\"\/', '\"' + oeis_url, s)
        if browse is None:
            if format == 'guess':
                if embedded():
                    return self.links(format='html')
                else:
                    return self.links(format='url')
            elif format == 'raw':
                return FancyTuple(self._fields['H'])
            elif format == 'html':
                return HtmlFragment(FancyTuple([url_absolute(_) for _ in self._fields['H']]))
            elif format == 'url':
                url_list = flatten([_urls(url_absolute(string)) for string in self._fields['H']])
                return FancyTuple(url_list)
        else:
            import webbrowser
            url_list = flatten([_urls(url_absolute(string)) for string in self._fields['H']])
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
            'F(n) = ((1+sqrt(5))^n-(1-sqrt(5))^n)/(2^n*sqrt(5)).'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.formulas()
            0: For n big enough, s(n+1) - s(n) = 0.

        """
        return FancyTuple(self._fields['F'])

    def cross_references(self, fetch=False):
        r"""
        Return a tuple of cross references associated to the sequence
        ``self``.

        INPUT:

        - ``fetch`` - boolean (default: ``False``).

        OUTPUT:

        - if ``fetch`` is ``False``, return a list of OEIS IDs (strings).
        - if ``fetch`` if ``True``, return a tuple of OEIS sequences.

        EXAMPLES::

            sage: nbalanced = oeis("A005598") ; nbalanced   # optional -- internet
            A005598: a(n)=1+sum((n-i+1)*phi(i),i=1..n).

            sage: nbalanced.cross_references()              # optional -- internet
            ('A049703', 'A049695', 'A103116', 'A000010')

            sage: nbalanced.cross_references(fetch=True)    # optional -- internet
            0: A049703: a(0) = 0; for n>0, a(n) = A005598(n)/2.
            1: A049695: Array T read by diagonals; T(i,j)=number of nonnegative slopes of lines determined by 2 lattice points in [ 0,i ] X [ 0,j ] if i>0; T(0,j)=1 if j>0; T(0,0)=0.
            2: A103116: A005598(n) - 1.
            3: A000010: Euler totient function phi(n): count numbers <= n and prime to n.

            sage: phi = _[3]                                # optional -- internet

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.cross_references()
            ('A000042', 'A000024')
        """
        ref_list = re.findall('A[0-9]{6}', " ".join(self._fields['Y']))
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
            A039834: a(n+2) = -a(n+1)+a(n) (signed Fibonacci numbers); or Fibonacci numbers (A000045) extended to negative indices.

            sage: sfibo.extensions_or_errors()[0]       # optional -- internet
            'Signs corrected by _Len Smiley_ and _N. J. A. Sloane_.'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.extensions_or_errors()
            0: This sequence does not contain errors.

        """
        return FancyTuple(self._fields['E'])

    def examples(self):
        r"""
        Return a tuple of examples associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: c = oeis(1203) ; c                    # optional -- internet
            A001203: Continued fraction expansion of Pi.

            sage: c.examples()                          # optional -- internet
            0: Pi = 3.1415926535897932384...
            1:    = 3 + 1/(7 + 1/(15 + 1/(1 + 1/(292 + ...))))
            2:    = [a_0; a_1, a_2, a_3, ...] = [3; 7, 15, 1, 292, ...]

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.examples()
            0: s(42) + s(43) = 0.
        """
        return FancyTuple(self._fields['e'])

    def comments(self):
        r"""
        Return a tuple of comments associated to the sequence ``self``.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.comments()[:3]                      # optional -- internet
            0: Also sometimes called Lamé's sequence.
            1: F(n+2) = number of binary sequences of length n that have no consecutive 0's.
            2: F(n+2) = number of subsets of {1,2,...,n} that contain no consecutive integers.

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.comments()
            0: 42 is the product of the first 4 prime numbers, except 5 and perhaps 1.
            1: Apart from that, i have no comment.
        """
        return FancyTuple(self._fields['C'])

    def url(self):
        r"""
        Return the URL of the page associated to the sequence ``self``.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: f = oeis(45) ; f                      # optional -- internet
            A000045: Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.

            sage: f.url()                               # optional -- internet
            'http://oeis.org/A000045'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.url()
            'http://oeis.org/A999999'
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
            FORMULAS
            ...
            OFFSETS
            (0, 1)
            <BLANKLINE>
            URL
            http://oeis.org/A012345
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
            The opposite of twice the characteristic sequence of 42 plus ...
            FIRST TERMS
            (1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
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
            if embedded() and s == 'links':
                print(re.sub('_', ' ', s).upper())
                getattr(self, s)()
                print('\n')
            else:
                result = getattr(self, s)()
                if result != '' and result != ('',) and result != ():
                    print(re.sub('_', ' ', s).upper())
                    print(str(result) + '\n')

    def programs(self, language='other'):
        r"""
        Returns programs implementing the sequence ``self`` in the given ``language``.

        INPUT:

        - ``language`` - string (default: 'other') - the language of the
          program. Current values are: 'maple', 'mathematica' and 'other'.

        OUTPUT:

        - tuple of strings (with fancy formatting).

        .. TODO:: ask OEIS to add a "Sage program" field in the database ;)

        EXAMPLES::

            sage: ee = oeis('A001113') ; ee             # optional -- internet
            A001113: Decimal expansion of e.

            sage: ee.programs()[0]                      # optional -- internet
            '(PARI) { default(realprecision, 50080); x=exp(1); for (n=1, 50000, d=floor(x); x=(x-d)*10; write("b001113.txt", n, " ", d)); } \\\\ _Harry J. Smith_, Apr 15 2009'

        TESTS::

            sage: s = oeis._imaginary_sequence()
            sage: s.programs()
            0: (Python)
            1: def A999999(n):
            2:     assert(isinstance(n, (int, Integer))), "n must be an integer."
            3:     if n < 38:
            4:         raise ValueError("The value %s is not accepted." %str(n)))
            5:     elif n == 42:
            6:         return -1
            7:     else:
            8:         return 1

            sage: s.programs('maple')
            0: Do not even try, Maple is not able to produce such a sequence.

            sage: s.programs('mathematica')
            0: Mathematica neither.
        """
        if language == "maple":
            return FancyTuple(self._fields['p'])
        elif language == "mathematica":
            return FancyTuple(self._fields['t'])
        else:
            return FancyTuple(self._fields['o'])


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
        Prints the tuple with one value per line, each line begins with the
        index of the value in ``self``.

        EXAMPLES::

            sage: from sage.databases.oeis import FancyTuple
            sage: t = FancyTuple(['zero', 'one', 'two', 'three', 4]) ; t
            0: zero
            1: one
            2: two
            3: three
            4: 4
        """
        length = len(str(len(self) - 1))
        return '\n'.join((('{0:>%d}' % length).format(str(i)) + ': ' + str(self[i]) for i in range(len(self))))

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
            sage: t
            ('\xc3\xa9', '\xc3\xa8', '\xc3\xa0', '\xc3\xa7')
            sage: FancyTuple(t)[2:4]
            0: à
            1: ç
        """
        return FancyTuple(tuple(self).__getslice__(i, j))

oeis = OEIS()
