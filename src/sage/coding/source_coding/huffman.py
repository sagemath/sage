r"""
Huffman Encoding

This module implements functionalities relating to Huffman encoding and
decoding.

AUTHOR:

- Nathann Cohen (2010-05): initial version.


Classes and functions
=====================
"""

###########################################################################
# Copyright (c) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# http://www.gnu.org/licenses/
###########################################################################
from __future__ import print_function

import six
from sage.structure.sage_object import SageObject

###########################################################################
#
# Helper functions
#
###########################################################################

def frequency_table(string):
    r"""
    Return the frequency table corresponding to the given string.

    INPUT:

    - ``string`` -- a string of symbols over some alphabet.

    OUTPUT:

    - A table of frequency of each unique symbol in ``string``. If ``string``
      is an empty string, return an empty table.

    EXAMPLES:

    The frequency table of a non-empty string::

        sage: from sage.coding.source_coding.huffman import frequency_table
        sage: str = "Stop counting my characters!"
        sage: T = sorted(frequency_table(str).items())
        sage: for symbol, code in T:
        ....:     print("{} {}".format(symbol, code))
          3
        ! 1
        S 1
        a 2
        c 3
        e 1
        g 1
        h 1
        i 1
        m 1
        n 2
        o 2
        p 1
        r 2
        s 1
        t 3
        u 1
        y 1

    The frequency of an empty string::

        sage: frequency_table("")
        {}
    """
    d = {}
    for s in string:
        d[s] = d.get(s, 0) + 1
    return d

class Huffman(SageObject):
    r"""
    This class implements the basic functionalities of Huffman codes.

    It can build a Huffman code from a given string, or from the information
    of a dictionary associating to each key (the elements of the alphabet) a
    weight (most of the time, a probability value or a number of occurrences).

    INPUT:

    - ``source`` -- can be either

        - A string from which the Huffman encoding should be created.

        - A dictionary that associates to each symbol of an alphabet a numeric
          value. If we consider the frequency of each alphabetic symbol, then
          ``source`` is considered as the frequency table of the alphabet with
          each numeric (non-negative integer) value being the number of
          occurrences of a symbol. The numeric values can also represent weights
          of the symbols. In that case, the numeric values are not necessarily
          integers, but can be real numbers.

    In order to construct a Huffman code for an alphabet, we use exactly one of
    the following methods:

    #. Let ``source`` be a string of symbols over an alphabet and feed
       ``source`` to the constructor of this class. Based on the input string, a
       frequency table is constructed that contains the frequency of each unique
       symbol in ``source``. The alphabet in question is then all the unique
       symbols in ``source``. A significant implication of this is that any
       subsequent string that we want to encode must contain only symbols that
       can be found in ``source``.

    #. Let ``source`` be the frequency table of an alphabet. We can feed this
       table to the constructor of this class. The table ``source`` can be a
       table of frequencies or a table of weights.

    Examples::

        sage: from sage.coding.source_coding.huffman import Huffman, frequency_table
        sage: h1 = Huffman("There once was a french fry")
        sage: for letter, code in h1.encoding_table().iteritems():
        ....:     print("'{}' : {}".format(letter, code))
        'a' : 0111
        ' ' : 00
        'c' : 1010
        'e' : 100
        'f' : 1011
        'h' : 1100
        'o' : 11100
        'n' : 1101
        's' : 11101
        'r' : 010
        'T' : 11110
        'w' : 11111
        'y' : 0110

    We can obtain the same result by "training" the Huffman code with the
    following table of frequency::

        sage: ft = frequency_table("There once was a french fry"); ft
        {' ': 5,
         'T': 1,
         'a': 2,
         'c': 2,
         'e': 4,
         'f': 2,
         'h': 2,
         'n': 2,
         'o': 1,
         'r': 3,
         's': 1,
         'w': 1,
         'y': 1}
        sage: h2 = Huffman(ft)

    Once ``h1`` has been trained, and hence possesses an encoding table,
    it is possible to obtain the Huffman encoding of any string
    (possibly the same) using this code::

        sage: encoded = h1.encode("There once was a french fry"); encoded
        '11110110010001010000111001101101010000111110111111010001110010110101001101101011000010110100110'

    We can decode the above encoded string in the following way::

        sage: h1.decode(encoded)
        'There once was a french fry'

    Obviously, if we try to decode a string using a Huffman instance which
    has been trained on a different sample (and hence has a different encoding
    table), we are likely to get some random-looking string::

        sage: h3 = Huffman("There once were two french fries")
        sage: h3.decode(encoded)
        ' wehnefetrhft ne ewrowrirTc'

    This does not look like our original string.

    Instead of using frequency, we can assign weights to each alphabetic
    symbol::

        sage: from sage.coding.source_coding.huffman import Huffman
        sage: T = {"a":45, "b":13, "c":12, "d":16, "e":9, "f":5}
        sage: H = Huffman(T)
        sage: L = ["deaf", "bead", "fab", "bee"]
        sage: E = []
        sage: for e in L:
        ....:     E.append(H.encode(e))
        ....:     print(E[-1])
        111110101100
        10111010111
        11000101
        10111011101
        sage: D = []
        sage: for e in E:
        ....:     D.append(H.decode(e))
        ....:     print(D[-1])
        deaf
        bead
        fab
        bee
        sage: D == L
        True
    """

    def __init__(self, source):
        r"""
        Constructor for Huffman.

        See the docstring of this class for full documentation.

        EXAMPLES::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)

        TESTS:

        Feeding anything else than a string or a dictionary::

            sage: Huffman(Graph())
            Traceback (most recent call last):
            ...
            ValueError: Input must be either a string or a dictionary.
        """

        # alphabetic symbol to Huffman encoding translation table
        self._character_to_code = []
        # Huffman binary tree
        self._tree = None
        # index of each alphabetic symbol
        self._index = None

        if isinstance(source, six.string_types):
            self._build_code(frequency_table(source))
        elif isinstance(source, dict):
            self._build_code(source)
        else:
            raise ValueError("Input must be either a string or a dictionary.")

    def _build_code_from_tree(self, tree, d, prefix):
        r"""
        Builds the Huffman code corresponding to a given tree and prefix.

        INPUT:

        - ``tree`` -- integer, or list of size `2`

        - ``d`` -- the dictionary to fill

        - ``prefix`` (string) -- binary string which is the prefix
          of any element of the tree

        EXAMPLES::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: d = {}
            sage: h._build_code_from_tree(h._tree, d, prefix="")
        """
        # This is really a recursive construction of a Huffman code. By
        # feeding this class a sufficiently large alphabet, it is possible to
        # exceed the maximum recursion depth and hence result in a RuntimeError.
        try:
            self._build_code_from_tree(tree[0],
                                       d,
                                       prefix="".join([prefix, "0"]))
            self._build_code_from_tree(tree[1],
                                       d,
                                       prefix="".join([prefix, "1"]))
        except TypeError:
            d[tree] = prefix

    def _build_code(self, dic):
        r"""
        Constructs a Huffman code corresponding to an alphabet with the given
        weight table.

        INPUT:

        - ``dic`` -- a dictionary that associates to each symbol of an alphabet
          a numeric value. If we consider the frequency of each alphabetic
          symbol, then ``dic`` is considered as the frequency table of the
          alphabet with each numeric (non-negative integer) value being the
          number of occurrences of a symbol. The numeric values can also
          represent weights of the symbols. In that case, the numeric values
          are not necessarily integers, but can be real numbers. In general,
          we refer to ``dic`` as a weight table.

        EXAMPLE::

            sage: from sage.coding.source_coding.huffman import Huffman, frequency_table
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: d = {}
            sage: h._build_code(frequency_table(str))
        """
        from heapq import heappush, heappop
        heap = []
        # Each alphabetic symbol is now represented by an element with
        # weight w and index i.
        for i, (s, w) in enumerate(dic.items()):
            heappush(heap, (w, i))
        for i in range(1, len(dic)):
            weight_a, node_a = heappop(heap)
            weight_b, node_b = heappop(heap)
            heappush(heap, (weight_a + weight_b, [node_a, node_b]))
        # dictionary of symbol to Huffman encoding
        d = {}
        self._tree = heap[0][1]
        # Build the binary tree of a Huffman code, where the root of the tree
        # is associated with the empty string.
        self._build_code_from_tree(self._tree, d, prefix="")
        self._index = dict((i, s) for i, (s, w) in enumerate(dic.items()))
        self._character_to_code = dict(
            (s, d[i]) for i, (s, w) in enumerate(dic.items()))

    def encode(self, string):
        r"""
        Encode the given string based on the current encoding table.

        INPUT:

        - ``string`` -- a string of symbols over an alphabet.

        OUTPUT:

        - A Huffman encoding of ``string``.

        EXAMPLES:

        This is how a string is encoded and then decoded::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: encoded = h.encode(str); encoded
            '00000110100010101011000011101010011100101010011011011100111101110010110100001011011111000001110101010001010110011010111111011001110100101000111110010011011100101011100000110001100101000101110101111101110110011000101011000111111101101111010010111001110100011'
            sage: h.decode(encoded)
            'Sage is my most favorite general purpose computer algebra system'
        """
        if self._character_to_code:
            return "".join((self._character_to_code[x] for x in string))

    def decode(self, string):
        r"""
        Decode the given string using the current encoding table.

        INPUT:

        - ``string`` -- a string of Huffman encodings.

        OUTPUT:

        - The Huffman decoding of ``string``.

        EXAMPLES:

        This is how a string is encoded and then decoded::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: encoded = h.encode(str); encoded
            '00000110100010101011000011101010011100101010011011011100111101110010110100001011011111000001110101010001010110011010111111011001110100101000111110010011011100101011100000110001100101000101110101111101110110011000101011000111111101101111010010111001110100011'
            sage: h.decode(encoded)
            'Sage is my most favorite general purpose computer algebra system'

        TESTS:

        Of course, the string one tries to decode has to be a binary one. If
        not, an exception is raised::

            sage: h.decode('I clearly am not a binary string')
            Traceback (most recent call last):
            ...
            ValueError: Input must be a binary string.
        """
        # This traverses the whole Huffman binary tree in order to work out
        # the symbol represented by a stream of binaries. This method of
        # decoding is really slow. A faster method is needed.
        # TODO: faster decoding implementation
        chars = []
        tree = self._tree
        index = self._index
        for i in string:
            if i == "0":
                tree = tree[0]
            elif i == "1":
                tree = tree[1]
            else:
                raise ValueError("Input must be a binary string.")
            if not isinstance(tree, list):
                chars.append(index[tree])
                tree = self._tree
        return "".join(chars)

    def encoding_table(self):
        r"""
        Returns the current encoding table.

        INPUT:

        - None.

        OUTPUT:

        - A dictionary associating an alphabetic symbol to a Huffman encoding.

        EXAMPLES::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: T = sorted(h.encoding_table().items())
            sage: for symbol, code in T:
            ....:     print("{} {}".format(symbol, code))
              101
            S 00000
            a 1101
            b 110001
            c 110000
            e 010
            f 110010
            g 0001
            i 10000
            l 10011
            m 0011
            n 110011
            o 0110
            p 0010
            r 1111
            s 1110
            t 0111
            u 10001
            v 00001
            y 10010
        """
        return self._character_to_code.copy()

    def tree(self):
        r"""
        Returns the Huffman tree corresponding to the current encoding.

        INPUT:

        - None.

        OUTPUT:

        - The binary tree representing a Huffman code.

        EXAMPLES::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: str = "Sage is my most favorite general purpose computer algebra system"
            sage: h = Huffman(str)
            sage: T = h.tree(); T
            Digraph on 39 vertices
            sage: T.show(figsize=[20,20])
            <BLANKLINE>
        """
        from sage.graphs.digraph import DiGraph
        g = DiGraph()
        g.add_edges(self._generate_edges(self._tree))
        return g

    def _generate_edges(self, tree, parent="", bit=""):
        """
        Generate the edges of the given Huffman tree.

        INPUT:

        - ``tree`` -- a Huffman binary tree.

        - ``parent`` -- (default: empty string) a parent vertex with exactly
          two children.

        - ``bit`` -- (default: empty string) the bit signifying either the
          left or right branch. The bit "0" denotes the left branch and "1"
          denotes the right branch.

        OUTPUT:

        - An edge list of the Huffman binary tree.

        EXAMPLES::

            sage: from sage.coding.source_coding.huffman import Huffman
            sage: H = Huffman("Sage")
            sage: T = H.tree()
            sage: T.edges(labels=None)
            [('0', 'S: 01'), ('0', 'a: 00'), ('1', 'e: 10'), ('1', 'g: 11'), ('root', '0'), ('root', '1')]
        """
        if parent == "":
            u = "root"
        else:
            u = parent
        s = "".join([parent, bit])
        try:
            left = self._generate_edges(tree[0], parent=s, bit="0")
            right = self._generate_edges(tree[1], parent=s, bit="1")
            L = [(u, s)] if s != "" else []
            return left + right + L
        except TypeError:
            return [(u, "".join([self.decode(s), ": ", s]))]
