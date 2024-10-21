# sage_setup: distribution = sagemath-combinat
# sage.doctest: needs sage.combinat sage.modules
r"""
Word paths

This module implements word paths, which is an application of Combinatorics
on Words to Discrete Geometry. A word path is the representation of a word
as a discrete path in a vector space using a one-to-one correspondence
between the alphabet and a set of vectors called steps. Many problems
surrounding 2d lattice polygons (such as questions of self-intersection,
area, inertia moment, etc.) can be solved in linear time (linear in the
length of the perimeter) using theory from Combinatorics on Words.

On the square grid, the encoding of a path using a four-letter alphabet
(for East, North, West and South directions) is also known as the Freeman
chain code [1,2] (see [3] for further reading).

AUTHORS:

- Arnaud Bergeron (2008) : Initial version, path on the square grid

- Sébastien Labbé (2009-01-14) : New classes and hierarchy, doc and functions.

EXAMPLES:

The combinatorial class of all paths defined over three given steps::

    sage: P = WordPaths('abc', steps=[(1,2), (-3,4), (0,-3)]); P
    Word Paths over 3 steps

This defines a one-to-one correspondence between alphabet and steps::

    sage: d = P.letters_to_steps()
    sage: sorted(d.items())
    [('a', (1, 2)), ('b', (-3, 4)), ('c', (0, -3))]

Creation of a path from the combinatorial class P defined above::

    sage: p = P('abaccba'); p
    Path: abaccba

Many functions can be used on p: the coordinates of its trajectory,
ask whether p is a closed path, plot it and many other::

    sage: list(p.points())
    [(0, 0), (1, 2), (-2, 6), (-1, 8), (-1, 5), (-1, 2), (-4, 6), (-3, 8)]
    sage: p.is_closed()
    False
    sage: p.plot()                                                                      # needs sage.plot
    Graphics object consisting of 3 graphics primitives

To obtain a list of all the available word path specific functions,
use ``help(p)``::

    sage: help(p)
    Help on FiniteWordPath_2d_str in module sage.combinat.words.paths object:
    ...
    Methods inherited from FiniteWordPath_2d:
    ...
    Methods inherited from FiniteWordPath_all:
    ...

Since p is a finite word, many functions from the word library are available::

    sage: p.crochemore_factorization()
    (a, b, a, c, c, ba)
    sage: p.is_palindrome()
    False
    sage: p[:3]
    Path: aba
    sage: len(p)
    7

P also herits many functions from Words::

    sage: P = WordPaths('rs', steps=[(1,2), (-1,4)]); P
    Word Paths over 2 steps
    sage: P.alphabet()
    {'r', 's'}
    sage: list(P.iterate_by_length(3))
    [Path: rrr,
     Path: rrs,
     Path: rsr,
     Path: rss,
     Path: srr,
     Path: srs,
     Path: ssr,
     Path: sss]

When the number of given steps is half the size of alphabet, the
opposite of vectors are used::

    sage: P = WordPaths('abcd', [(1,0), (0,1)])
    sage: sorted(P.letters_to_steps().items())
    [('a', (1, 0)), ('b', (0, 1)), ('c', (-1, 0)), ('d', (0, -1))]

Some built-in combinatorial classes of paths::

    sage: P = WordPaths('abAB', steps='square_grid'); P
    Word Paths on the square grid

::

    sage: D = WordPaths('()', steps='dyck'); D
    Finite Dyck paths
    sage: d = D('()()()(())'); d
    Path: ()()()(())
    sage: d.plot()                                                                      # needs sage.plot
    Graphics object consisting of 3 graphics primitives

::

    sage: P = WordPaths('abcdef', steps='triangle_grid')
    sage: p = P('babaddefadabcadefaadfafabacdefa')
    sage: p.plot()                                                                      # needs sage.plot
    Graphics object consisting of 3 graphics primitives

Vector steps may be in more than 2 dimensions::

    sage: d = [(1,0,0), (0,1,0), (0,0,1)]
    sage: P = WordPaths(alphabet='abc', steps=d); P
    Word Paths over 3 steps
    sage: p = P('abcabcabcabcaabacabcababcacbabacacabcaccbcac')
    sage: p.plot()                                                                      # needs sage.plot
    Graphics3d Object

::

    sage: d = [(1,3,5,1), (-5,1,-6,0), (0,0,1,9), (4,2,-1,0)]
    sage: P = WordPaths(alphabet='rstu', steps=d); P
    Word Paths over 4 steps
    sage: p = P('rtusuusususuturrsust'); p
    Path: rtusuusususuturrsust
    sage: p.end_point()
    (5, 31, -26, 30)

::

    sage: CubePaths = WordPaths('abcABC', steps='cube_grid'); CubePaths
    Word Paths on the cube grid
    sage: CubePaths('abcabaabcabAAAAA').plot()                                          # needs sage.plot
    Graphics3d Object

The input data may be a str, a list, a tuple,
a callable or a finite iterator::

    sage: P = WordPaths([0, 1, 2, 3])
    sage: P([0,1,2,3,2,1,2,3,2])
    Path: 012321232
    sage: P((0,1,2,3,2,1,2,3,2))
    Path: 012321232
    sage: P(lambda n:n%4, length=10)
    Path: 0123012301
    sage: P(iter([0,3,2,1]), length='finite')
    Path: 0321

REFERENCES:

- [1] Freeman, H.: *On the encoding of arbitrary geometric configurations*.
  IRE Trans. Electronic Computer 10 (1961) 260-268.
- [2] Freeman, H.: *Boundary encoding and processing*. In Lipkin, B., Rosenfeld,
  A., eds.: Picture Processing and Psychopictorics, Academic Press, New York
  (1970) 241-266.
- [3] Braquelaire, J.P., Vialard, A.: *Euclidean paths: A new representation of
  boundary of discrete regions*. Graphical Models and Image Processing 61 (1999)
  16-43.
- [4] :wikipedia:`Regular_tiling`
- [5] :wikipedia:`Dyck_word`
"""
# ****************************************************************************
#       Copyright (C) 2008 Arnaud bergeron <abergeron@gmail.coms>,
#       Copyright (C) 2009 Sebastien Labbe <slabqc@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.words.words import FiniteWords
from sage.combinat.words.word import FiniteWord_class
from sage.combinat.words.alphabet import build_alphabet
from sage.misc.lazy_import import lazy_import
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from .word_datatypes import (WordDatatype_str,
                             WordDatatype_list,
                             WordDatatype_tuple)
# WordDatatype_cpp_basic_string)

from .word_infinite_datatypes import (
    WordDatatype_iter_with_caching,
    WordDatatype_iter,
    WordDatatype_callable_with_caching,
    WordDatatype_callable)
from sage.matrix.constructor import vector_on_axis_rotation_matrix

lazy_import("sage.plot.all", ["arrow", "line", "polygon", "point", "Graphics"])
lazy_import('sage.rings.number_field.number_field', 'QuadraticField')


#######################################################################
#                                                                     #
#                         WordPaths function                          #
#                                                                     #
#######################################################################

def WordPaths(alphabet, steps=None):
    r"""
    Return the combinatorial class of paths of the given type of steps.

    INPUT:

    - ``alphabet`` -- ordered alphabet

    - ``steps`` -- (default: ``None``) it can be one of the following:

      - an iterable ordered container of as many vectors as there are
        letters in the alphabet. The vectors are associated to the letters
        according to their order in steps. The vectors can be a tuple or
        anything that can be passed to vector function.

      - an iterable ordered container of k vectors where k is half the
        size of alphabet. The vectors and their opposites are associated
        to the letters according to their order in steps (given vectors
        first, opposite vectors after).

      - ``None`` -- in this case, the type of steps are guessed from the
        length of alphabet

      - ``'square_grid'`` or ``'square'`` -- (default when size of alphabet is 4)
        The order is : East, North, West, South.

      - ``'triangle_grid'`` or ``'triangle'``

      - ``'hexagonal_grid'`` or ``'hexagon'`` -- (default when size of alphabet is 6)

      - ``'cube_grid'`` or ``'cube'``

      - ``'north_east'``, ``'ne'`` or ``'NE'`` -- (the default when size of alphabet is 2)

      - ``'dyck'``

    OUTPUT: the combinatorial class of all paths of the given type

    EXAMPLES:

    The steps can be given explicitly::

        sage: WordPaths('abc', steps=[(1,2), (-1,4), (0,-3)])
        Word Paths over 3 steps

    Different type of input alphabet::

        sage: WordPaths(range(3), steps=[(1,2), (-1,4), (0,-3)])
        Word Paths over 3 steps
        sage: WordPaths(['cric','crac','croc'], steps=[(1,2), (1,4), (0,3)])
        Word Paths over 3 steps

    Directions can be in three dimensions as well::

        sage: WordPaths('ab', steps=[(1,2,2),(-1,4,2)])
        Word Paths over 2 steps

    When the number of given steps is half the size of alphabet, the
    opposite of vectors are used::

        sage: P = WordPaths('abcd', [(1,0), (0,1)])
        sage: P
        Word Paths over 4 steps
        sage: sorted(P.letters_to_steps().items())
        [('a', (1, 0)), ('b', (0, 1)), ('c', (-1, 0)), ('d', (0, -1))]

    When no steps are given, default classes are returned::

        sage: WordPaths('ab')
        Word Paths in North and East steps
        sage: WordPaths(range(4))
        Word Paths on the square grid
        sage: WordPaths(range(6))
        Word Paths on the hexagonal grid

    There are many type of built-in steps...

    On a two letters alphabet::

        sage: WordPaths('ab', steps='north_east')
        Word Paths in North and East steps
        sage: WordPaths('()', steps='dyck')
        Finite Dyck paths

    On a four letters alphabet::

        sage: WordPaths('ruld', steps='square_grid')
        Word Paths on the square grid

    On a six letters alphabet::

        sage: WordPaths('abcdef', steps='hexagonal_grid')
        Word Paths on the hexagonal grid
        sage: WordPaths('abcdef', steps='triangle_grid')
        Word Paths on the triangle grid
        sage: WordPaths('abcdef', steps='cube_grid')
        Word Paths on the cube grid

    TESTS::

        sage: WordPaths(range(5))
        Traceback (most recent call last):
        ...
        TypeError: Unable to make a class WordPaths from {0, 1, 2, 3, 4}
        sage: WordPaths('abAB', steps='square_gridd')
        Traceback (most recent call last):
        ...
        TypeError: Unknown type of steps : square_gridd
    """
    # Construction of the alphabet
    alphabet = build_alphabet(alphabet)

    # If no steps are given, they are guessed from the alphabet
    if steps is None:
        if alphabet.cardinality() == 2:
            steps = 'north_east'
        elif alphabet.cardinality() == 4:
            steps = 'square_grid'
        elif alphabet.cardinality() == 6:
            steps = 'hexagonal_grid'
        else:
            raise TypeError("Unable to make a class WordPaths from %s" % alphabet)

    # Return the class of WordPaths according to the given type of paths
    if isinstance(steps, str):
        if steps in ('square_grid', 'square'):
            return WordPaths_square_grid(alphabet=alphabet)
        elif steps in ('triangle_grid', 'triangle'):
            return WordPaths_triangle_grid(alphabet=alphabet)
        elif steps in ('hexagonal_grid', 'hexagon'):
            return WordPaths_hexagonal_grid(alphabet=alphabet)
        elif steps in ('cube_grid', 'cube'):
            return WordPaths_cube_grid(alphabet=alphabet)
        elif steps in ('north_east', 'ne', 'NE'):
            return WordPaths_north_east(alphabet=alphabet)
        elif steps == 'dyck':
            return WordPaths_dyck(alphabet=alphabet)
        else:
            raise TypeError("Unknown type of steps : %s" % steps)
    else:
        return WordPaths_all(alphabet=alphabet, steps=steps)


#######################################################################
#                                                                     #
#                  Combinatorial classes of word paths                #
#                                                                     #
#######################################################################

class WordPaths_all(FiniteWords):
    r"""
    The combinatorial class of all paths, i.e of all words over
    an alphabet where each letter is mapped to a step (a vector).
    """
    def __init__(self, alphabet, steps):
        r"""
        INPUT:

        - ``alphabet`` -- an ordered alphabet

        - ``steps`` -- an iterable (of same length as alphabet or half the
          length of alphabet) of ordered vectors

        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_all
            sage: d = ((1,1), (-1,1), (1,-1), (-1,-1))
            sage: P = WordPaths_all('abAB', d); P
            Word Paths over 4 steps
            sage: P == loads(dumps(P))
            True

        If size of alphabet is twice the number of steps, then opposite
        vectors are used for the second part of the alphabet::

            sage: WordPaths('abcd',[(2,1),(2,4)])
            Word Paths over 4 steps
            sage: _.letters_to_steps()
            {'a': (2, 1), 'b': (2, 4), 'c': (-2, -1), 'd': (-2, -4)}

        TESTS::

            sage: from sage.combinat.words.paths import WordPaths_all
            sage: d = ((1,1), (-1,1), (1,-1), (-1,-1))
            sage: WordPaths_all('abA', d)
            Traceback (most recent call last):
            ...
            TypeError: size of steps (=4) must equal the size of alphabet (=3) or half the size of alphabet

            sage: d = ((1,1), 1)
            sage: WordPaths_all('ab', d)
            Traceback (most recent call last):
            ...
            ValueError: cannot make vectors from steps

            sage: d = ((1,1), (-1,1,0))
            sage: WordPaths_all('ab', d)
            Traceback (most recent call last):
            ...
            ValueError: cannot make summable vectors from steps
        """
        # Construction of the words class
        FiniteWords.__init__(self, alphabet)
        alphabet = self.alphabet()

        # Checking the size of alphabet and steps
        ls = len(steps)
        la = alphabet.cardinality()
        if la != ls and la != 2 * ls:
            raise TypeError("size of steps (=%s) must equal the size of alphabet (=%s) or half the size of alphabet" % (len(steps), alphabet.cardinality()))

        # Construction of the steps
        from sage.structure.element import Vector
        if all(isinstance(x, Vector) for x in steps):
            vsteps = steps
        else:
            try:
                vsteps = [vector(s) for s in steps]
            except (TypeError):
                raise ValueError("cannot make vectors from steps")
        try:
            s = sum(vsteps)
        except (TypeError, AttributeError):
            raise ValueError("cannot make summable vectors from steps")

        # Complete vsteps with the opposite vectors if needed
        if la == 2 * ls:
            vsteps += [-v for v in vsteps]

        self._steps = dict(zip(alphabet, vsteps))
        self._vector_space = s.parent()

    def __eq__(self, other):
        r"""
        TESTS::

            sage: W1 = WordPaths(['a','b'], [vector((0,1)), vector((0,2))])
            sage: W2 = WordPaths(['a','b'], [vector((0,1)), vector((0,2))])
            sage: W3 = WordPaths(['a','b'], [vector((0,2)), vector((1,0))])
            sage: W1 == W2
            True
            sage: W1 == W3
            False
        """
        return self is other or (type(self) is type(other) and
                                 self.alphabet() == other.alphabet() and
                                 self.vector_space() == other.vector_space() and
                                 self.letters_to_steps() == other.letters_to_steps())

    def __ne__(self, other):
        r"""
        TESTS::

            sage: W1 = WordPaths(['a','b'], [vector((0,1)), vector((0,2))])
            sage: W2 = WordPaths(['a','b'], [vector((0,1)), vector((0,2))])
            sage: W3 = WordPaths(['a','b'], [vector((0,2)), vector((1,0))])
            sage: W1 != W2
            False
            sage: W1 != W3
            True
        """
        return not (self == other)

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the elements of ``self``.

        The word may be finite (infinite or of unknown length is not supported
        yet).
        Its data may be str, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.
        The dimension of the path may be 1, 2, 3 or more.

        TESTS::

            sage: d = WordPaths('ab',steps=[(1,2),(3,4)])._element_classes
            sage: type(d)
            <class 'dict'>
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_2d_tuple'>

        ::

            sage: d = WordPaths('ab',steps=[(1,2,3),(3,4,5)])._element_classes
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_3d_tuple'>

        ::

            sage: steps = [(1,2,3,4),(3,4,5,6)]
            sage: d = WordPaths('ab',steps=steps)._element_classes
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_all_tuple'>

        ::

            sage: d = WordPaths('ab',steps=[(1,),(3,)])._element_classes
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_all_tuple'>
        """
        dimension = self._vector_space.dimension()
        if dimension == 2:
            return {
                'list': FiniteWordPath_2d_list,
                'str': FiniteWordPath_2d_str,
                'tuple': FiniteWordPath_2d_tuple,
                'callable_with_caching': FiniteWordPath_2d_callable_with_caching,
                'callable': FiniteWordPath_2d_callable,
                'iter_with_caching': FiniteWordPath_2d_iter_with_caching,
                'iter': FiniteWordPath_2d_iter,
            }
        elif dimension == 3:
            return {
                'list': FiniteWordPath_3d_list,
                'str': FiniteWordPath_3d_str,
                'tuple': FiniteWordPath_3d_tuple,
                'callable_with_caching': FiniteWordPath_3d_callable_with_caching,
                'callable': FiniteWordPath_3d_callable,
                'iter_with_caching': FiniteWordPath_3d_iter_with_caching,
                'iter': FiniteWordPath_3d_iter,
            }
        else:
            return {
                'list': FiniteWordPath_all_list,
                'str': FiniteWordPath_all_str,
                'tuple': FiniteWordPath_all_tuple,
                'callable_with_caching': FiniteWordPath_all_callable_with_caching,
                'callable': FiniteWordPath_all_callable,
                'iter_with_caching': FiniteWordPath_all_iter_with_caching,
                'iter': FiniteWordPath_all_iter,
            }

    def __repr__(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_all
            sage: d = (vector((1,1)), vector((-1,1)), vector((1,-1)), vector((-1,-1)))
            sage: WordPaths_all('abAB',d).__repr__()
            'Word Paths over 4 steps'
        """
        return "Word Paths over %s steps" % self.alphabet().cardinality()

    def letters_to_steps(self) -> dict:
        r"""
        Return the dictionary mapping letters to vectors (steps).

        EXAMPLES::

            sage: d = WordPaths('ab').letters_to_steps()
            sage: sorted(d.items())
            [('a', (0, 1)), ('b', (1, 0))]
            sage: d = WordPaths('abcd').letters_to_steps()
            sage: sorted(d.items())
            [('a', (1, 0)), ('b', (0, 1)), ('c', (-1, 0)), ('d', (0, -1))]
            sage: d = WordPaths('abcdef').letters_to_steps()
            sage: sorted(d.items())
            [('a', (1, 0)),
             ('b', (1/2, 1/2*sqrt3)),
             ('c', (-1/2, 1/2*sqrt3)),
             ('d', (-1, 0)),
             ('e', (-1/2, -1/2*sqrt3)),
             ('f', (1/2, -1/2*sqrt3))]
        """
        return self._steps

    def vector_space(self):
        r"""
        Return the vector space over which the steps of the paths are defined.

        EXAMPLES::

            sage: WordPaths('ab',steps='dyck').vector_space()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: WordPaths('ab',steps='north_east').vector_space()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: WordPaths('abcd',steps='square_grid').vector_space()
            Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: WordPaths('abcdef',steps='hexagonal_grid').vector_space()
            Vector space of dimension 2 over Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?
            sage: WordPaths('abcdef',steps='cube_grid').vector_space()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
            sage: WordPaths('abcdef',steps='triangle_grid').vector_space()
            Vector space of dimension 2 over Number Field in sqrt3 with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?
        """
        return self._vector_space


class WordPaths_square_grid(WordPaths_all):
    r"""
    The combinatorial class of all paths on the square grid.
    """
    def __init__(self, alphabet):
        r"""
        The combinatorial class of all finite paths on the square grid.

        INPUT:

        - ``alphabet`` -- ordered alphabet of length 4; the order for the steps
          is : East, North, West, South

        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_square_grid
            sage: P = WordPaths_square_grid('abAB'); P
            Word Paths on the square grid
            sage: P == loads(dumps(P))
            True
        """
        # Construction of the steps
        d = [(1, 0), (0, 1), (-1, 0), (0, -1)]

        # Construction of the class
        super().__init__(alphabet, steps=d)

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the elements of ``self``.

        The word may be finite (infinite or of unknown length is not supported
        yet).
        Its data may be str, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        TESTS::

            sage: d = WordPaths('abcd')._element_classes
            sage: type(d)
            <class 'dict'>
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_square_grid_tuple'>
        """
        return {
            'list': FiniteWordPath_square_grid_list,
            'str': FiniteWordPath_square_grid_str,
            'tuple': FiniteWordPath_square_grid_tuple,
            'callable_with_caching': FiniteWordPath_square_grid_callable_with_caching,
            'callable': FiniteWordPath_square_grid_callable,
            'iter_with_caching': FiniteWordPath_square_grid_iter_with_caching,
            'iter': FiniteWordPath_square_grid_iter,
        }

    def __repr__(self) -> str:
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_square_grid
            sage: WordPaths_square_grid('abAB').__repr__()
            'Word Paths on the square grid'
        """
        return "Word Paths on the square grid"


class WordPaths_triangle_grid(WordPaths_all):
    r"""
    The combinatorial class of all paths on the triangle grid.
    """
    def __init__(self, alphabet):
        r"""
        The combinatorial class of all finite paths on the triangle grid.

        INPUT:

        - ``alphabet`` -- ordered alphabet of length 6. The order for the steps
          is : Right, Up-Right, Up-Left, Left, Down-Left, Down-Right.

        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_triangle_grid
            sage: P = WordPaths_triangle_grid('abcdef'); P
            Word Paths on the triangle grid
            sage: P == loads(dumps(P))
            True
        """
        K = QuadraticField(3, 'sqrt3')
        sqrt3 = K.gen()

        # Construction of the steps
        d = (vector(K, (1, 0)),
             vector(K, (ZZ(1) / ZZ(2), sqrt3 / 2)),
             vector(K, (ZZ(-1) / ZZ(2), sqrt3 / 2)),
             vector(K, (-1, 0)),
             vector(K, (ZZ(-1) / ZZ(2), -sqrt3 / 2)),
             vector(K, (ZZ(1) / ZZ(2), -sqrt3 / 2)))

        # Construction of the class
        super().__init__(alphabet, steps=d)

        self._infinite_word_class = None
        self._finite_word_class = FiniteWordPath_triangle_grid

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the elements of ``self``.

        The word may be finite (infinite or of unknown length is not supported
        yet).
        Its data may be str, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        TESTS::

            sage: d = WordPaths('abcdef', steps='triangle')._element_classes
            sage: len(d)
            7
            sage: type(d)
            <class 'dict'>
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_triangle_grid_tuple'>
        """
        return {
            'list': FiniteWordPath_triangle_grid_list,
            'str': FiniteWordPath_triangle_grid_str,
            'tuple': FiniteWordPath_triangle_grid_tuple,
            'callable_with_caching': FiniteWordPath_triangle_grid_callable_with_caching,
            'callable': FiniteWordPath_triangle_grid_callable,
            'iter_with_caching': FiniteWordPath_triangle_grid_iter_with_caching,
            'iter': FiniteWordPath_triangle_grid_iter,
        }

    def __repr__(self) -> str:
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_triangle_grid
            sage: WordPaths_triangle_grid('abcdef').__repr__()
            'Word Paths on the triangle grid'
        """
        return "Word Paths on the triangle grid"


class WordPaths_hexagonal_grid(WordPaths_triangle_grid):
    r"""
    The combinatorial class of all paths on the hexagonal grid.
    """
    def __init__(self, alphabet):
        r"""
        The combinatorial class of all finite paths on the hexagonal grid.

        INPUT:

        - ``alphabet`` -- ordered alphabet of length 6. The order for the steps
          is : Right, Up-Right, Up-Left, Left, Down-Left, Down-Right.

        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_hexagonal_grid
            sage: P = WordPaths_hexagonal_grid('abcdef'); P
            Word Paths on the hexagonal grid
            sage: P == loads(dumps(P))
            True
        """
        # Construction of the class
        super().__init__(alphabet)

        self._infinite_word_class = None
        self._finite_word_class = FiniteWordPath_hexagonal_grid

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the elements of ``self``.

        The word may be finite (infinite or of unknown length is not supported
        yet).
        Its data may be str, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        TESTS::

            sage: d = WordPaths('abcdef', steps='hexagon')._element_classes
            sage: type(d)
            <class 'dict'>
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_hexagonal_grid_tuple'>
        """
        return {
            'list': FiniteWordPath_hexagonal_grid_list,
            'str': FiniteWordPath_hexagonal_grid_str,
            'tuple': FiniteWordPath_hexagonal_grid_tuple,
            'callable_with_caching': FiniteWordPath_hexagonal_grid_callable_with_caching,
            'callable': FiniteWordPath_hexagonal_grid_callable,
            'iter_with_caching': FiniteWordPath_hexagonal_grid_iter_with_caching,
            'iter': FiniteWordPath_hexagonal_grid_iter,
        }

    def __repr__(self) -> str:
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_hexagonal_grid
            sage: WordPaths_hexagonal_grid('abcdef').__repr__()
            'Word Paths on the hexagonal grid'
        """
        return "Word Paths on the hexagonal grid"


class WordPaths_cube_grid(WordPaths_all):
    r"""
    The combinatorial class of all paths on the cube grid.
    """
    def __init__(self, alphabet):
        r"""
        The combinatorial class of all finite paths on the cube grid.

        INPUT:

        - ``alphabet`` -- ordered alphabet of length 6. The order for
          the steps is `e_x, e_y, e_z, -e_x, -e_y, -e_z`, where `e_v`
          denotes the canonical basis.

        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_cube_grid
            sage: P = WordPaths_cube_grid('abcABC'); P
            Word Paths on the cube grid
            sage: P == loads(dumps(P))
            True
        """
        # Construction of the class
        d = [(1, 0, 0), (0, 1, 0), (0, 0, 1),
             (-1, 0, 0), (0, -1, 0), (0, 0, -1)]
        super().__init__(alphabet, steps=d)
        self._infinite_word_class = None
        self._finite_word_class = FiniteWordPath_cube_grid

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the elements of ``self``.

        The word may be finite (infinite or of unknown length is not supported
        yet).
        Its data may be str, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        TESTS::

            sage: d = WordPaths('abcdef', steps='cube')._element_classes
            sage: type(d)
            <class 'dict'>
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_cube_grid_tuple'>
        """
        return {'list': FiniteWordPath_cube_grid_list,
                'str': FiniteWordPath_cube_grid_str,
                'tuple': FiniteWordPath_cube_grid_tuple,
                'callable_with_caching': FiniteWordPath_cube_grid_callable_with_caching,
                'callable': FiniteWordPath_cube_grid_callable,
                'iter_with_caching': FiniteWordPath_cube_grid_iter_with_caching,
                'iter': FiniteWordPath_cube_grid_iter,
                }

    def __repr__(self) -> str:
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_cube_grid
            sage: WordPaths_cube_grid('abcABC').__repr__()
            'Word Paths on the cube grid'
        """
        return "Word Paths on the cube grid"


class WordPaths_dyck(WordPaths_all):
    r"""
    The combinatorial class of all Dyck paths.
    """
    def __init__(self, alphabet):
        r"""
        The combinatorial class of all finite Dyck paths.

        INPUT:

        - ``alphabet`` -- ordered alphabet of length 2. The order for the steps
          is : (1,1), (1,-1)

        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_dyck
            sage: P = WordPaths_dyck('[]'); P
            Finite Dyck paths
            sage: P == loads(dumps(P))
            True
        """
        # Construction of the class
        d = [(1, 1), (1, -1)]
        super().__init__(alphabet, steps=d)

        self._infinite_word_class = None
        self._finite_word_class = FiniteWordPath_dyck

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the elements of ``self``.

        The word may be finite (infinite or of unknown length is not supported
        yet).
        Its data may be str, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        TESTS::

            sage: d = WordPaths('ab', steps='dyck')._element_classes
            sage: type(d)
            <class 'dict'>
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_dyck_tuple'>
        """
        return {'list': FiniteWordPath_dyck_list,
                'str': FiniteWordPath_dyck_str,
                'tuple': FiniteWordPath_dyck_tuple,
                'callable_with_caching': FiniteWordPath_dyck_callable_with_caching,
                'callable': FiniteWordPath_dyck_callable,
                'iter_with_caching': FiniteWordPath_dyck_iter_with_caching,
                'iter': FiniteWordPath_dyck_iter,
                }

    def __repr__(self) -> str:
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_dyck
            sage: WordPaths_dyck('()').__repr__()
            'Finite Dyck paths'
        """
        return "Finite Dyck paths"


class WordPaths_north_east(WordPaths_all):
    r"""
    The combinatorial class of all paths using North and East directions.
    """
    def __init__(self, alphabet):
        r"""
        The combinatorial class of all finite paths using only north and east
        steps on the square grid.

        INPUT:

        - ``alphabet`` -- ordered alphabet of length 2. The order for the steps
          is North, East

        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_north_east
            sage: P = WordPaths_north_east('ab'); P
            Word Paths in North and East steps
            sage: P == loads(dumps(P))
            True
        """
        # Construction of the class
        d = [(0, 1), (1, 0)]
        super().__init__(alphabet, steps=d)
        self._infinite_word_class = None
        self._finite_word_class = FiniteWordPath_north_east

    @lazy_attribute
    def _element_classes(self):
        r"""
        Return a dictionary that gives the class of the elements of ``self``.

        The word may be finite (infinite or of unknown length is not supported
        yet).
        Its data may be str, list, tuple, a callable or an iterable.
        For callable and iterable, the data may be cached.

        TESTS::

            sage: d = WordPaths('ab', steps='NE')._element_classes
            sage: type(d)
            <class 'dict'>
            sage: len(d)
            7
            sage: d['tuple']
            <class 'sage.combinat.words.paths.FiniteWordPath_north_east_tuple'>
        """
        return {'list': FiniteWordPath_north_east_list,
                'str': FiniteWordPath_north_east_str,
                'tuple': FiniteWordPath_north_east_tuple,
                'callable_with_caching': FiniteWordPath_north_east_callable_with_caching,
                'callable': FiniteWordPath_north_east_callable,
                'iter_with_caching': FiniteWordPath_north_east_iter_with_caching,
                'iter': FiniteWordPath_north_east_iter,
                }

    def __repr__(self) -> str:
        r"""
        EXAMPLES::

            sage: from sage.combinat.words.paths import WordPaths_north_east
            sage: WordPaths_north_east('ab').__repr__()
            'Word Paths in North and East steps'
        """
        return "Word Paths in North and East steps"


#######################################################################
#                                                                     #
#                    Abstract word path classes                       #
#                        (all, 2d, 3d, ...)                           #
#                                                                     #
#######################################################################

class FiniteWordPath_all(SageObject):
    def _repr_(self) -> str:
        r"""
        Return a string representation of this path.

        EXAMPLES::

            sage: F = WordPaths('ab',[(1,0,0,0),(0,1,0,0)]); F
            Word Paths over 2 steps
            sage: f = F('ababab')
            sage: f._repr_()
            'Path: ababab'
        """
        return "Path: %s" % self.string_rep()

    def points(self, include_last=True):
        r"""
        Return an iterator yielding a list of points used to draw the path
        represented by this word.

        INPUT:

        - ``include_last`` -- boolean (default: ``True``); whether to include the
          last point

        EXAMPLES:

        A simple closed square::

            sage: P = WordPaths('abAB')
            sage: list(P('abAB').points())
            [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]

        A simple closed square without the last point::

            sage: list(P('abAB').points(include_last=False))
            [(0, 0), (1, 0), (1, 1), (0, 1)]

        ::

            sage: list(P('abaB').points())
            [(0, 0), (1, 0), (1, 1), (2, 1), (2, 0)]
        """
        curpt = self.start_point()
        yield curpt
        end = len(self) if include_last else -1
        for l in self[:end]:
            curpt += self.parent().letters_to_steps()[l]
            yield curpt

    def start_point(self):
        r"""
        Return the starting point of ``self``.

        OUTPUT: vector

        EXAMPLES::

            sage: WordPaths('abcdef')('abcdef').start_point()
            (0, 0)
            sage: WordPaths('abcdef', steps='cube_grid')('abcdef').start_point()
            (0, 0, 0)
            sage: P = WordPaths('ab', steps=[(1,0,0,0),(0,1,0,0)])
            sage: P('abbba').start_point()
            (0, 0, 0, 0)
        """
        return self.parent().vector_space()(0)

    @cached_method
    def end_point(self):
        r"""
        Return the end point of the path.

        EXAMPLES::

            sage: WordPaths('abcdef')('abababab').end_point()
            (6, 2*sqrt3)
            sage: WordPaths('abAB')('abababab').end_point()
            (4, 4)
            sage: P = WordPaths('abcABC', steps='cube_grid')
            sage: P('ababababCC').end_point()
            (4, 4, -2)
            sage: WordPaths('abcdef')('abcdef').end_point()
            (0, 0)
            sage: P = WordPaths('abc', steps=[(1,3,7,9),(-4,1,0,0),(0,32,1,8)])
            sage: P('abcabababacaacccbbcac').end_point()
            (-16, 254, 63, 128)
        """
        last = None
        for pt in self.points():
            last = pt
        return last

    def directive_vector(self):
        r"""
        Return the directive vector of ``self``.

        The directive vector is the vector starting at the start point
        and ending at the end point of the path ``self``.

        EXAMPLES::

            sage: WordPaths('abcdef')('abababab').directive_vector()
            (6, 2*sqrt3)
            sage: WordPaths('abAB')('abababab').directive_vector()
            (4, 4)
            sage: P = WordPaths('abcABC', steps='cube_grid')
            sage: P('ababababCC').directive_vector()
            (4, 4, -2)
            sage: WordPaths('abcdef')('abcdef').directive_vector()
            (0, 0)
            sage: P = WordPaths('abc', steps=[(1,3,7,9),(-4,1,0,0),(0,32,1,8)])
            sage: P('abcabababacaacccbbcac').directive_vector()
            (-16, 254, 63, 128)
        """
        return self.end_point() - self.start_point()

    def is_closed(self) -> bool:
        r"""
        Return ``True`` if the path is closed.

        A path is closed if the origin and the end of
        the path are equal.

        EXAMPLES::

            sage: P = WordPaths('abcd', steps=[(1,0),(0,1),(-1,0),(0,-1)])
            sage: P('abcd').is_closed()
            True
            sage: P('abc').is_closed()
            False
            sage: P().is_closed()
            True
            sage: P('aacacc').is_closed()
            True
        """
        return self.start_point() == self.end_point()

    def is_simple(self) -> bool:
        r"""
        Return ``True`` if the path is simple.

        A path is simple if all its points are
        distinct.

        If the path is closed, the last point is not considered.

        EXAMPLES::

            sage: P = WordPaths('abcdef',steps='triangle_grid');P
            Word Paths on the triangle grid
            sage: P('abc').is_simple()
            True
            sage: P('abcde').is_simple()
            True
            sage: P('abcdef').is_simple()
            True
            sage: P('ad').is_simple()
            True
            sage: P('aabdee').is_simple()
            False
        """
        n = 0
        s = set()
        include_last = not self.is_closed()
        for p in self.points(include_last=include_last):
            # We need the elements to have a common parent,
            # so we convert the points to immutable vectors.
            v = vector(p)
            v.set_immutable()
            s.add(v)
            n += 1
            if len(s) != n:
                return False
        return True

    def tikz_trajectory(self) -> str:
        r"""
        Return the trajectory of ``self`` as a ``tikz`` string.

        EXAMPLES::

            sage: P = WordPaths('abcdef')
            sage: p = P('abcde')
            sage: p.tikz_trajectory()
            '(0.000, 0.000) -- (1.00, 0.000) -- (1.50, 0.866) -- (1.00, 1.73) -- (0.000, 1.73) -- (-0.500, 0.866)'
        """
        from sage.misc.functional import N as n
        l = (str(tuple(n(x, digits=3) for x in pt)) for pt in self.points())
        return ' -- '.join(l)

    def projected_point_iterator(self, v=None, ring=None):
        r"""
        Return an iterator of the projection of the orbit points of the
        path into the space orthogonal to the given vector.

        INPUT:

        - ``v`` -- vector (default: ``None``); if ``None``, the directive
          vector (i.e. the end point minus starting point) of the path is
          considered

        - ``ring`` -- ring (default: ``None``); where to do the
          computations. If ``None``, RealField(53) is used.

        OUTPUT: iterator of points

        EXAMPLES:

        Projected points of the Rauzy fractal::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: D = s.fixed_point('1')
            sage: v = s.pisot_eigenvector_right()
            sage: P = WordPaths('123',[(1,0,0),(0,1,0),(0,0,1)])
            sage: w = P(D[:200])
            sage: it = w.projected_point_iterator(v)
            sage: for i in range(6): next(it)
            (0.000000000000000, 0.000000000000000)
            (-0.526233343362516, 0.000000000000000)
            (0.220830337618112, -0.477656250512816)
            (-0.305403005744404, -0.477656250512816)
            (0.100767309386062, 0.400890564600664)
            (-0.425466033976454, 0.400890564600664)

        Projected points of a 2d path::

            sage: P = WordPaths('ab','ne')
            sage: p = P('aabbabbab')
            sage: it = p.projected_point_iterator(ring=RealField(20))
            sage: for i in range(8): next(it)
            (0.00000)
            (0.78087)
            (1.5617)
            (0.93704)
            (0.31235)
            (1.0932)
            (0.46852)
            (-0.15617)
        """
        if v is None:
            v = self.directive_vector()
        if ring is None:
            ring = RR
        R = vector_on_axis_rotation_matrix(v, 0, ring=ring)[1:]
        for q in self.points():
            yield R * q

    def plot_projection(self, v=None, letters=None, color=None, ring=None,
                        size=12, kind='right'):
        r"""
        Return an image of the projection of the successive points of the
        path into the space orthogonal to the given vector.

        INPUT:

        - ``self`` -- a word path in a 3 or 4 dimension vector space

        - ``v`` -- vector (default: ``None``); if ``None``, the directive
          vector (i.e. the end point minus starting point) of the path is
          considered.

        - ``letters`` -- iterable (default: ``None``); of the letters
          to be projected. If ``None``, then all the letters are considered.

        - ``color`` -- dictionary (default: ``None``); of the letters
          mapped to colors. If ``None``, automatic colors are chosen.

        - ``ring`` -- ring (default: ``None``); where to do the
          computations. If ``None``, RealField(53) is used.

        - ``size`` -- number (default: ``12``); size of the points

        - ``kind`` -- string (default: ``'right'``); either
          ``'right'`` or ``'left'``. The color of a letter is given to the
          projected prefix to the right or the left of the letter.

        OUTPUT: 2d or 3d Graphic object

        EXAMPLES:

        The Rauzy fractal::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: D = s.fixed_point('1')
            sage: v = s.pisot_eigenvector_right()
            sage: P = WordPaths('123',[(1,0,0),(0,1,0),(0,0,1)])
            sage: w = P(D[:200])
            sage: w.plot_projection(v)  # long time (2s)
            Graphics object consisting of 200 graphics primitives

        In this case, the abelianized vector doesn't give a good
        projection::

            sage: w.plot_projection()  # long time (2s)
            Graphics object consisting of 200 graphics primitives

        You can project only the letters you want::

            sage: w.plot_projection(v, letters='12')  # long time (2s)
            Graphics object consisting of 168 graphics primitives

        You can increase or decrease the precision of the computations by
        changing the ring of the projection matrix::

            sage: w.plot_projection(v, ring=RealField(20))  # long time (2s)
            Graphics object consisting of 200 graphics primitives

        You can change the size of the points::

            sage: w.plot_projection(v, size=30)  # long time (2s)
            Graphics object consisting of 200 graphics primitives

        You can assign the color of a letter to the projected prefix to the
        right or the left of the letter::

            sage: w.plot_projection(v, kind='left')  # long time (2s)
            Graphics object consisting of 200 graphics primitives

        To remove the axis, do like this::

            sage: r = w.plot_projection(v)                                              # needs sage.plot
            sage: r.axes(False)                                                         # needs sage.plot
            sage: r                             # long time (2s)                        # needs sage.plot
            Graphics object consisting of 200 graphics primitives

        You can assign different colors to each letter::

            sage: color = {'1': 'purple', '2': (.2,.3,.4), '3': 'magenta'}
            sage: w.plot_projection(v, color=color)     # long time (2s)                # needs sage.plot
            Graphics object consisting of 200 graphics primitives

        The 3d-Rauzy fractal::

            sage: s = WordMorphism('1->12,2->13,3->14,4->1')
            sage: D = s.fixed_point('1')
            sage: v = s.pisot_eigenvector_right()
            sage: P = WordPaths('1234',[(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)])
            sage: w = P(D[:200])
            sage: w.plot_projection(v)                                                  # needs sage.plot
            Graphics3d Object

        The dimension of vector space of the parent must be 3 or 4::

            sage: P = WordPaths('ab', [(1, 0), (0, 1)])
            sage: p = P('aabbabbab')
            sage: p.plot_projection()                                                   # needs sage.plot
            Traceback (most recent call last):
            ...
            TypeError: The dimension of the vector space (=2) must be 3 or 4
        """
        dimension = self.parent().vector_space().dimension()
        if dimension not in (3, 4):
            msg = "The dimension of the vector space (=%s) must be 3 or 4" % dimension
            raise TypeError(msg)
        if letters is None:
            letters = self.parent().alphabet()
        if color is None:
            from sage.plot.all import hue
            A = self.parent().alphabet()
            color = {a: hue(A.rank(a) / float(A.cardinality())) for a in A}
        it = self.projected_point_iterator(v, ring=ring)
        if kind == 'right':
            next(it)
        elif kind != 'left':
            raise ValueError('unknown value for kind (=%s)' % kind)
        tout = [point([c], color=color[a], size=size)
                for a, c in zip(self, it) if a in letters]
        return sum(tout)

    def projected_path(self, v=None, ring=None):
        r"""
        Return the path projected into the space orthogonal to the given
        vector.

        INPUT:

        - ``v`` -- vector (default: ``None``); if ``None``, the directive
          vector (i.e. the end point minus starting point) of the path is
          considered.

        - ``ring`` -- ring (default: ``None``); where to do the
          computations. If ``None``, RealField(53) is used.

        OUTPUT: word path

        EXAMPLES:

        The projected path of the tribonacci word::

            sage: s = WordMorphism('1->12,2->13,3->1')
            sage: D = s.fixed_point('1')
            sage: v = s.pisot_eigenvector_right()
            sage: P = WordPaths('123',[(1,0,0),(0,1,0),(0,0,1)])
            sage: w = P(D[:1000])
            sage: p = w.projected_path(v)
            sage: p
            Path: 1213121121312121312112131213121121312121...
            sage: p[:20].plot()                                                         # needs sage.plot
            Graphics object consisting of 3 graphics primitives

        The ``ring`` argument allows to change the precision of the
        projected steps::

            sage: p = w.projected_path(v, RealField(10))
            sage: p
            Path: 1213121121312121312112131213121121312121...
            sage: p.parent().letters_to_steps()
            {'1': (-0.53, 0.00), '2': (0.75, -0.48), '3': (0.41, 0.88)}
        """
        if v is None:
            v = self.directive_vector()
        if ring is None:
            ring = RR
        R = vector_on_axis_rotation_matrix(v, 0, ring=ring)[1:]
        d = self.parent().letters_to_steps()
        A = self.parent().alphabet()
        nvvectors = [R * d[a] for a in A]
        projected_parent = WordPaths(A, nvvectors)
        return projected_parent(self)

    def is_tangent(self):
        r"""
        The is_tangent() method, which is implemented for words, has
        an extended meaning for word paths, which is not implemented yet.

        TESTS::

            sage: WordPaths('ab')('abbab').is_tangent()
            Traceback (most recent call last):
            ...
            NotImplementedError

        AUTHOR:

        -   Thierry Monteil
        """
        raise NotImplementedError


class FiniteWordPath_2d(FiniteWordPath_all):
    def plot(self, pathoptions={"rgbcolor": 'red', "thickness": 3},
             fill=True, filloptions={"rgbcolor": 'red', "alpha": 0.2},
             startpoint=True, startoptions={"rgbcolor": 'red', "pointsize": 100},
             endarrow=True, arrowoptions={"rgbcolor": 'red', "arrowsize": 20, "width": 3},
             gridlines=False, gridoptions={}):
        r"""
        Return a 2d Graphics illustrating the path.

        INPUT:

        - ``pathoptions`` -- (dict,
          default:dict(rgbcolor='red',thickness=3)), options for the
          path drawing

        - ``fill`` -- boolean (default: ``True``); if fill is ``True`` and if
          the path is closed, the inside is colored

        - ``filloptions`` -- (dict,
          default:dict(rgbcolor='red',alpha=0.2)), options for the
          inside filling

        - ``startpoint`` -- boolean (default: ``True``); draw the start point?

        - ``startoptions`` -- (dict,
          default:dict(rgbcolor='red',pointsize=100)) options for the
          start point drawing

        - ``endarrow`` -- boolean (default: ``True``); draw an arrow end at the end?

        - ``arrowoptions`` -- (dict,
          default:dict(rgbcolor='red',arrowsize=20, width=3)) options
          for the end point arrow

        - ``gridlines`` -- boolean (default: ``False``); show gridlines?

        - ``gridoptions`` -- (dict, default: {}), options for the gridlines

        EXAMPLES:

        A non closed path on the square grid::

            sage: P = WordPaths('abAB')
            sage: P('abababAABAB').plot()                                               # needs sage.plot
            Graphics object consisting of 3 graphics primitives

        A closed path on the square grid::

            sage: P('abababAABABB').plot()                                              # needs sage.plot
            Graphics object consisting of 4 graphics primitives

        A Dyck path::

            sage: P = WordPaths('()', steps='dyck')
            sage: P('()()()((()))').plot()                                              # needs sage.plot
            Graphics object consisting of 3 graphics primitives

        A path in the triangle grid::

            sage: P = WordPaths('abcdef', steps='triangle_grid')
            sage: P('abcdedededefab').plot()                                            # needs sage.plot
            Graphics object consisting of 3 graphics primitives

        A polygon of length 220 that tiles the plane in two ways::

            sage: P = WordPaths('abAB')
            sage: P('aBababAbabaBaBABaBabaBaBABAbABABaBabaBaBABaBababAbabaBaBABaBabaBaBABAbABABaBABAbAbabAbABABaBABAbABABaBabaBaBABAbABABaBABAbAbabAbABAbAbabaBababAbABAbAbabAbABABaBABAbAbabAbABAbAbabaBababAbabaBaBABaBababAbabaBababAbABAbAbab').plot()  # needs sage.plot
            Graphics object consisting of 4 graphics primitives

        With gridlines::

            sage: P('ababababab').plot(gridlines=True)                                  # needs sage.plot

        TESTS::

            sage: P = WordPaths('abAB')
            sage: P().plot()                                                            # needs sage.plot
            Graphics object consisting of 3 graphics primitives
            sage: sum(map(plot,map(P,['a','A','b','B'])))                               # needs sage.plot
            Graphics object consisting of 12 graphics primitives
        """
        G = Graphics()
        pts = list(self.points())

        ####################
        ####################
        # FIXME Bug: plot needs float for coordinates
        ####################
        ####################
        pts = [[RR(i) for i in x] for x in pts]

        # Inside
        if fill and self.is_closed():
            G += polygon(pts, **filloptions)

        # Startpoint
        if startpoint:
            G += point(pts[0], **startoptions)

        # The path itself
        if endarrow and not self.is_empty():
            G += line(pts[:-1], **pathoptions)
            G += arrow(pts[-2], pts[-1], **arrowoptions)
        else:
            G += line(pts, **pathoptions)

        G.axes(False)
        G.set_aspect_ratio(1)

        # gridlines
        # ############## BUG ##############
        # Gridlines doesn't work fine.
        # It should be gridlines="integers"
        # ############## BUG ##############
        if gridlines:
            G = G.show(gridlines=True, **gridoptions)

        return G

    def animate(self):
        r"""
        Return an animation object illustrating the path growing step by step.

        EXAMPLES::

            sage: P = WordPaths('abAB')
            sage: p = P('aaababbb')
            sage: a = p.animate(); print(a)                                             # needs sage.plot
            Animation with 9 frames
            sage: show(a)                       # long time, optional - imagemagick, needs sage.plot
            sage: show(a, delay=35, iterations=3)       # long time, optional - imagemagick, needs sage.plot

        ::

            sage: P = WordPaths('abcdef',steps='triangle')
            sage: p =  P('abcdef')
            sage: a = p.animate(); print(a)                                             # needs sage.plot
            Animation with 8 frames
            sage: show(a)                       # long time, optional - imagemagick, needs sage.plot

        If the path is closed, the plain polygon is added at the end of the
        animation::

            sage: P = WordPaths('abAB')
            sage: p = P('ababAbABABaB')
            sage: a = p.animate(); print(a)                                             # needs sage.plot
            Animation with 14 frames
            sage: show(a)                       # long time, optional - imagemagick, needs sage.plot

        Another example illustrating a Fibonacci tile::

            sage: w = words.fibonacci_tile(2)
            sage: a = w.animate(); print(a)                                             # needs sage.plot
            Animation with 54 frames
            sage: show(a)                       # long time, optional - imagemagick, needs sage.plot

        The first 4 Fibonacci tiles in an animation::

            sage: # needs sage.plot
            sage: a = words.fibonacci_tile(0).animate()
            sage: b = words.fibonacci_tile(1).animate()
            sage: c = words.fibonacci_tile(2).animate()
            sage: d = words.fibonacci_tile(3).animate()
            sage: print(a*b*c*d)
            Animation with 296 frames
            sage: show(a*b*c*d)                 # long time, optional - imagemagick

        .. NOTE::

            If ImageMagick is not installed, you will get an error
            message like this::

               convert: not found

               Error: ImageMagick does not appear to be installed. Saving an
               animation to a GIF file or displaying an animation requires
               ImageMagick, so please install it and try again.

            See www.imagemagick.org, for example.
        """
        from sage.plot.all import line, polygon, animate

        pts = list(self.points())

        ####################
        ####################
        # Bug: plot needs float for coordinates
        ####################
        ####################
        pts = [[RR(i) for i in x] for x in pts]

        images = [line(pts[:i]) for i in range(1, len(pts) + 1)]

        if self.is_closed():
            images.append(polygon(pts))

        # Get the window of the last image
        last_image = images[-1]
        kwds = {}
        kwds['xmin'] = last_image.xmin()
        kwds['xmax'] = last_image.xmax()
        kwds['ymin'] = last_image.ymin()
        kwds['ymax'] = last_image.ymax()
        kwds['aspect_ratio'] = 1
        kwds['axes'] = False

        return animate(images, **kwds)

    def plot_directive_vector(self, options={"rgbcolor": 'blue'}):
        r"""
        Return an arrow 2d graphics that goes from the start of the path
        to the end.

        INPUT:

        - ``options`` -- dictionary, default: {'rgbcolor': 'blue'} graphic
          options for the arrow

        If the start is the same as the end, a single point is returned.

        EXAMPLES::

            sage: P = WordPaths('abcd'); P
            Word Paths on the square grid
            sage: p = P('aaaccaccacacacaccccccbbdd'); p
            Path: aaaccaccacacacaccccccbbdd
            sage: R = p.plot() + p.plot_directive_vector()                              # needs sage.plot
            sage: R.axes(False)                                                         # needs sage.plot
            sage: R.set_aspect_ratio(1)                                                 # needs sage.plot
            sage: R.plot()                                                              # needs sage.plot
            Graphics object consisting of 4 graphics primitives

        TESTS:

        A closed path::

            sage: P('acbd').plot_directive_vector()                                     # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        start = self.start_point()
        end = self.end_point()
        if (start == end):
            G = point(start, pointsize=10, **options)
        else:
            G = arrow(start, end, **options)
        G.axes(False)
        G.set_aspect_ratio(1)
        return G

    def area(self):
        r"""
        Return the area of a closed path.

        INPUT:

        - ``self`` -- a closed path

        EXAMPLES::

            sage: P = WordPaths('abcd',steps=[(1,1),(-1,1),(-1,-1),(1,-1)])
            sage: p = P('abcd')
            sage: p.area()          #todo: not implemented
            2
        """
        if not self.is_closed():
            raise TypeError("the path must be closed to compute its area")
        return NotImplemented

    def height(self):
        r"""
        Return the height of ``self``.

        The height of a `2d`-path is merely the difference
        between the highest and the lowest `y`-coordinate of each
        points traced by it.

        OUTPUT: nonnegative real number

        EXAMPLES::

            sage: Freeman = WordPaths('abAB')
            sage: Freeman('aababaabbbAA').height()
            5

        The function is well-defined if ``self`` is not simple or close::

            sage: Freeman('aabAAB').height()
            1
            sage: Freeman('abbABa').height()
            2

        This works for any `2d`-paths::

            sage: Paths = WordPaths('ab', steps=[(1,0),(1,1)])
            sage: p = Paths('abbaa')
            sage: p.height()
            2
            sage: DyckPaths = WordPaths('ab', steps='dyck')
            sage: p = DyckPaths('abaabb')
            sage: p.height()
            2
            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.height()
            2.59807621135332
        """
        return self.ymax() - self.ymin()

    def height_vector(self):
        r"""
        Return the height at each point.

        EXAMPLES::

            sage: Paths = WordPaths('ab', steps=[(1,0),(0,1)])
            sage: p = Paths('abbba')
            sage: p.height_vector()
            [0, 0, 1, 2, 3, 3]
        """
        h_vec = []
        y_min = None
        y_max = None
        for _, y in self.points():
            if y_min is None:
                y_min = y
                y_max = y
            else:
                y_max = max(y, y_max)
                y_min = min(y, y_min)
            h_vec.append(y_max - y_min)
        return h_vec

    def width(self):
        r"""
        Return the width of ``self``.

        The height of a `2d`-path is merely the difference
        between the rightmost and the leftmost `x`-coordinate of each
        points traced by it.

        OUTPUT: nonnegative real number

        EXAMPLES::

            sage: Freeman = WordPaths('abAB')
            sage: Freeman('aababaabbbAA').width()
            5

        The function is well-defined if ``self`` is not simple or close::

            sage: Freeman('aabAAB').width()
            2
            sage: Freeman('abbABa').width()
            1

        This works for any `2d`-paths::

            sage: Paths = WordPaths('ab', steps=[(1,0),(1,1)])
            sage: p = Paths('abbaa')
            sage: p.width()
            5
            sage: DyckPaths = WordPaths('ab', steps='dyck')
            sage: p = DyckPaths('abaabb')
            sage: p.width()
            6
            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.width()
            4.50000000000000
        """
        return self.xmax() - self.xmin()

    def width_vector(self):
        r"""
        Return the width at each point.

        EXAMPLES::

            sage: Paths = WordPaths('ab', steps=[(1,0),(0,1)])
            sage: p = Paths('abbba')
            sage: p.width_vector()
            [0, 1, 1, 1, 1, 2]
        """
        w_vec = []
        x_min = None
        x_max = None
        for x, _ in self.points():
            if x_min is None:
                x_min = x
                x_max = x
            else:
                x_max = max(x, x_max)
                x_min = min(x, x_min)
            w_vec.append(x_max - x_min)
        return w_vec

    def xmin(self):
        r"""
        Return the minimum of the x-coordinates of the path.

        EXAMPLES::

            sage: P = WordPaths('0123')
            sage: p = P('0101013332')
            sage: p.xmin()
            0

        This works for any `2d`-paths::

            sage: Paths = WordPaths('ab', steps=[(1,0),(-1,1)])
            sage: p = Paths('abbba')
            sage: p.xmin()
            -2
            sage: DyckPaths = WordPaths('ab', steps='dyck')
            sage: p = DyckPaths('abaabb')
            sage: p.xmin()
            0
            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.xmin()
            0.000000000000000
        """
        return min(x for (x, _) in self.points())

    def ymin(self):
        r"""
        Return the minimum of the y-coordinates of the path.

        EXAMPLES::

            sage: P = WordPaths('0123')
            sage: p = P('0101013332')
            sage: p.ymin()
            0

        This works for any `2d`-paths::

            sage: Paths = WordPaths('ab', steps=[(1,-1),(-1,1)])
            sage: p = Paths('ababa')
            sage: p.ymin()
            -1
            sage: DyckPaths = WordPaths('ab', steps='dyck')
            sage: p = DyckPaths('abaabb')
            sage: p.ymin()
            0
            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.ymin()
            0.000000000000000
        """
        return min(y for (_, y) in self.points())

    def xmax(self):
        r"""
        Return the maximum of the x-coordinates of the path.

        EXAMPLES::

            sage: P = WordPaths('0123')
            sage: p = P('0101013332')
            sage: p.xmax()
            3

        This works for any `2d`-paths::

            sage: Paths = WordPaths('ab', steps=[(1,-1),(-1,1)])
            sage: p = Paths('ababa')
            sage: p.xmax()
            1
            sage: DyckPaths = WordPaths('ab', steps='dyck')
            sage: p = DyckPaths('abaabb')
            sage: p.xmax()
            6
            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.xmax()
            4.50000000000000
        """
        return max(x for (x, _) in self.points())

    def ymax(self):
        r"""
        Return the maximum of the y-coordinates of the path.

        EXAMPLES::

            sage: P = WordPaths('0123')
            sage: p = P('0101013332')
            sage: p.ymax()
            3

        This works for any `2d`-paths::

            sage: Paths = WordPaths('ab', steps=[(1,-1),(-1,1)])
            sage: p = Paths('ababa')
            sage: p.ymax()
            0
            sage: DyckPaths = WordPaths('ab', steps='dyck')
            sage: p = DyckPaths('abaabb')
            sage: p.ymax()
            2
            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.ymax()
            2.59807621135332
        """
        return max(y for (_, y) in self.points())


class FiniteWordPath_3d(FiniteWordPath_all):
    def plot(self, pathoptions={"rgbcolor": 'red', "arrow_head": True, "thickness": 3},
             startpoint=True, startoptions={"rgbcolor": 'red', "size": 10}):
        r"""
        INPUT:

        - ``pathoptions`` -- (dict, default:dict(rgbcolor='red',arrow_head=True,
          thickness=3)), options for the path drawing

        - ``startpoint`` -- boolean (default: ``True``); draw the start point?

        - ``startoptions`` -- (dict, default:dict(rgbcolor='red',size=10))
          options for the start point drawing

        EXAMPLES::

            sage: d = ( vector((1,3,2)), vector((2,-4,5)) )
            sage: P = WordPaths(alphabet='ab', steps=d); P
            Word Paths over 2 steps
            sage: p = P('ababab'); p
            Path: ababab
            sage: p.plot()                                                              # needs sage.plot
            Graphics3d Object

            sage: P = WordPaths('abcABC', steps='cube_grid')
            sage: p = P('abcabcAABBC')
            sage: p.plot()                                                              # needs sage.plot
            Graphics3d Object
        """
        # The following line seems not to work for 3d
        # G = Graphics()
        # so we draw to start a small almost invisible point instead:
        G = point([self.start_point()], size=1)
        pts = list(self.points())
        if startpoint:
            G += point([pts[0]], **startoptions)
        G += line(pts, **pathoptions)
        return G


#######################################################################
#                                                                     #
#                    Abstract word path classes                       #
#                (square grid, hexagonal grid, etc.)                  #
#                                                                     #
#######################################################################

class FiniteWordPath_square_grid(FiniteWordPath_2d):
    def is_closed(self) -> bool:
        r"""
        Return whether ``self`` represents a closed path.

        EXAMPLES::

            sage: P = WordPaths('abAB', steps='square_grid')
            sage: P('aA').is_closed()
            True
            sage: P('abAB').is_closed()
            True
            sage: P('ababAABB').is_closed()
            True
            sage: P('aaabbbAABB').is_closed()
            False
            sage: P('ab').is_closed()
            False
        """
        tab = self.abelian_vector()
        return tab[0] == tab[2] and tab[1] == tab[3]

    def area(self):
        r"""
        Return the area of a closed path.

        INPUT:

        - ``self`` -- a closed path

        EXAMPLES::

            sage: P = WordPaths('abAB', steps='square_grid')
            sage: P('abAB').area()
            1
            sage: P('aabbAABB').area()
            4
            sage: P('aabbABAB').area()
            3

        The area of the Fibonacci tiles::

            sage: [words.fibonacci_tile(i).area() for i in range(6)]
            [1, 5, 29, 169, 985, 5741]
            sage: [words.dual_fibonacci_tile(i).area() for i in range(6)]
            [1, 5, 29, 169, 985, 5741]
            sage: oeis(_)[0]                            # optional -- internet
            A001653: Numbers k such that 2*k^2 - 1 is a square.
            sage: _.first_terms()                       # optional -- internet
            (1,
             5,
             29,
             169,
             985,
             5741,
             33461,
             195025,
             1136689,
             6625109,
             38613965,
             225058681,
             1311738121,
             7645370045,
             44560482149,
             259717522849,
             1513744654945,
             8822750406821,
             51422757785981,
             299713796309065,
             1746860020068409,
             10181446324101389,
             59341817924539925)

        TESTS::

            sage: P = WordPaths('abAB', steps='square_grid')
            sage: P('a').area()
            Traceback (most recent call last):
            ...
            TypeError: the path must be closed to compute its area
        """
        if not self.is_closed():
            raise TypeError("the path must be closed to compute its area")
        return abs(self._area_vh())

    def _area_vh(self, x=0, y=0):
        r"""
        Return the area of ``self``, with starting point (x,y).

        This is using VH algorithm.

        INPUT:

        - x, y -- starting point (default: (0, 0))

        EXAMPLES::

            sage: P = WordPaths('abAB', steps='square_grid')
            sage: P('abAB')._area_vh()
            -1
            sage: P('aabbAABB')._area_vh()
            -4
            sage: P('aabbABAB')._area_vh()
            -3

        REFERENCES:

        Annie Lacasse Memoire.
        """
        area = 0
        a, b, A, B = self.parent().alphabet()

        for move in self:
            if move == b:
                area -= x
                y += 1
            elif move == B:
                area += x
                y -= 1
            elif move == a:
                area += y
                x += 1
            elif move == A:
                area -= y
                x -= 1
        return area // 2

    def is_simple(self) -> bool:
        r"""
        Return whether the path is simple.

        A path is simple if all its points are
        distinct.

        If the path is closed, the last point is not considered.

        .. NOTE::

            The linear algorithm described in the thesis of Xavier Provençal
            should be implemented here.

        EXAMPLES::

            sage: P = WordPaths('abAB', steps='square_grid')
            sage: P('abab').is_simple()
            True
            sage: P('abAB').is_simple()
            True
            sage: P('abA').is_simple()
            True
            sage: P('aabABB').is_simple()
            False
            sage: P().is_simple()
            True
            sage: P('A').is_simple()
            True
            sage: P('aA').is_simple()
            True
            sage: P('aaA').is_simple()
            False

        REFERENCES:

        - Provençal, X., *Combinatoires des mots, géometrie discrète et
          pavages*, Thèse de doctorat en Mathématiques, Montréal, UQAM,
          septembre 2008, 115 pages.
        """
        return super().is_simple()

    def tikz_trajectory(self) -> str:
        r"""
        Return the trajectory of ``self`` as a ``tikz`` string.

        EXAMPLES::

            sage: f = words.fibonacci_tile(1)
            sage: f.tikz_trajectory()
            '(0, 0) -- (0, -1) -- (-1, -1) -- (-1, -2) -- (0, -2) -- (0, -3) -- (1, -3) -- (1, -2) -- (2, -2) -- (2, -1) -- (1, -1) -- (1, 0) -- (0, 0)'
        """
        return ' -- '.join(map(str, self.points()))


class FiniteWordPath_triangle_grid(FiniteWordPath_2d):
    # Triangle grid paths are implemented with quadratic fields,
    # and the ordering of such elements is currently problematic:
    #
    #     sage: Q.<sqrt3> = QuadraticField(3)
    #     sage: sqrt3 > 0
    #     True
    #     sage: 0 < sqrt3
    #     False
    #     sage: max(2*sqrt3, sqrt3/10)
    #     1/10*sqrt3
    #
    # Therefore, the functions xmin(), xmax(), ymin() and ymax() are
    # redefined here with conversion to RR in order to avoid this problem
    def xmin(self):
        r"""
        Return the minimum of the x-coordinates of the path.

        EXAMPLES::

            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.xmin()
            0.000000000000000
            sage: w = WordPaths('abcABC', steps='triangle')('ABAcacacababababcbcbAC')
            sage: w.xmin()
            -3.00000000000000
        """
        return min(RR(x) for (x, _) in self.points())

    def ymin(self):
        r"""
        Return the minimum of the y-coordinates of the path.

        EXAMPLES::

            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.ymin()
            0.000000000000000
            sage: w = WordPaths('abcABC', steps='triangle')('ABAcacacababababcbcbAC')
            sage: w.ymin()
            -0.866025403784439
        """
        return min(RR(y) for (_, y) in self.points())

    def xmax(self):
        r"""
        Return the maximum of the x-coordinates of the path.

        EXAMPLES::

            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.xmax()
            4.50000000000000
            sage: w = WordPaths('abcABC', steps='triangle')('ABAcacacababababcbcbAC')
            sage: w.xmax()
            4.00000000000000
        """
        return max(RR(x) for (x, _) in self.points())

    def ymax(self):
        r"""
        Return the maximum of the y-coordinates of the path.

        EXAMPLES::

            sage: w = WordPaths('abcABC', steps='triangle')('ababcaaBC')
            sage: w.ymax()
            2.59807621135332
            sage: w = WordPaths('abcABC', steps='triangle')('ABAcacacababababcbcbAC')
            sage: w.ymax()
            8.66025403784439
        """
        return max(RR(y) for (_, y) in self.points())


# TODO: faire une verification du mot pour etre sur hexagonal grid
class FiniteWordPath_hexagonal_grid(FiniteWordPath_triangle_grid):
    def __init__(self, parent, *args, **kwds):
        r"""
        INPUT:

        - ``parent`` -- a parent object inheriting from Words_all
          that has the alphabet attribute defined

        - ``*args``, ``**kwds`` -- arguments accepted by AbstractWord

        EXAMPLES::

            sage: F = WordPaths('abcdef', steps='hexagon'); F
            Word Paths on the hexagonal grid
            sage: f = F('aaabbbccddef'); f
            Path: aaabbbccddef

        ::

            sage: f == loads(dumps(f))
            True
        """
        super().__init__(parent, *args, **kwds)


class FiniteWordPath_cube_grid(FiniteWordPath_3d):
    pass


class FiniteWordPath_north_east(FiniteWordPath_2d):
    pass


class FiniteWordPath_dyck(FiniteWordPath_2d):
    pass


#######################################################################
#                                                                     #
#                    Concrete word path classes                       #
#                                                                     #
#          It would be nice if those were created inline...           #
#       We must ask if Nicolas Thiery was able to convince Sage       #
#       people about this.                                            #
#                                                                     #
#######################################################################

# #### Finite paths ####

class FiniteWordPath_all_list(WordDatatype_list, FiniteWordPath_all, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths(['a','b'],[(1,2,0,0),(3,4,0,0)])
        sage: p = P(['a','b','a']);p
        Path: aba
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_all_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_all_str(WordDatatype_str, FiniteWordPath_all, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('ab',[(1,2,0,0),(3,4,0,0)])
        sage: p = P('aabbb'); p
        Path: aabbb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_all_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_all_tuple(WordDatatype_tuple, FiniteWordPath_all, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('ab',[(1,2,0,0),(3,4,0,0)])
        sage: p = P( ('a','b','b') ); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_all_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_all_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_all, FiniteWord_class):
    pass


class FiniteWordPath_all_iter(WordDatatype_iter, FiniteWordPath_all, FiniteWord_class):
    pass


class FiniteWordPath_all_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_all, FiniteWord_class):
    pass


class FiniteWordPath_all_callable(WordDatatype_callable, FiniteWordPath_all, FiniteWord_class):
    pass


# #### Finite paths on 2d ####

class FiniteWordPath_2d_list(WordDatatype_list, FiniteWordPath_2d, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths(['a','b'],[(1,2),(3,4)])
        sage: p = P(['a','b','a']);p
        Path: aba
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_2d_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_2d_str(WordDatatype_str, FiniteWordPath_2d, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths(['a','b'],[(1,2),(3,4)])
        sage: p = P('aba'); p
        Path: aba
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_2d_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_2d_tuple(WordDatatype_tuple, FiniteWordPath_2d, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths(['a','b'],[(1,2),(3,4)])
        sage: p = P(('a','b','a'));p
        Path: aba
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_2d_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_2d_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_2d, FiniteWord_class):
    pass


class FiniteWordPath_2d_iter(WordDatatype_iter, FiniteWordPath_2d, FiniteWord_class):
    pass


class FiniteWordPath_2d_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_2d, FiniteWord_class):
    pass


class FiniteWordPath_2d_callable(WordDatatype_callable, FiniteWordPath_2d, FiniteWord_class):
    pass


# #### Finite paths on 3d ####

class FiniteWordPath_3d_list(WordDatatype_list, FiniteWordPath_3d, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths(['a','b'],[(1,2,0),(3,4,0)])
        sage: p = P(['a','b','a']);p
        Path: aba
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_3d_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_3d_str(WordDatatype_str, FiniteWordPath_3d, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths(['a','b'],[(1,2,0),(3,4,0)])
        sage: p = P('aba'); p
        Path: aba
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_3d_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_3d_tuple(WordDatatype_tuple, FiniteWordPath_3d, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths(['a','b'],[(1,2,0),(3,4,0)])
        sage: p = P(('a','b','a'));p
        Path: aba
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_3d_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_3d_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_3d, FiniteWord_class):
    pass


class FiniteWordPath_3d_iter(WordDatatype_iter, FiniteWordPath_3d, FiniteWord_class):
    pass


class FiniteWordPath_3d_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_3d, FiniteWord_class):
    pass


class FiniteWordPath_3d_callable(WordDatatype_callable, FiniteWordPath_3d, FiniteWord_class):
    pass


# #### Finite paths on square grid ####

class FiniteWordPath_square_grid_list(WordDatatype_list, FiniteWordPath_square_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcd', steps='square')
        sage: p = P(['a','b','b']); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_square_grid_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_square_grid_str(WordDatatype_str, FiniteWordPath_square_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcd', steps='square')
        sage: p = P('abccc'); p
        Path: abccc
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_square_grid_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_square_grid_tuple(WordDatatype_tuple, FiniteWordPath_square_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcd', steps='square')
        sage: p = P(('a','b','b')); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_square_grid_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_square_grid_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_square_grid, FiniteWord_class):
    pass


class FiniteWordPath_square_grid_iter(WordDatatype_iter, FiniteWordPath_square_grid, FiniteWord_class):
    pass


class FiniteWordPath_square_grid_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_square_grid, FiniteWord_class):
    pass


class FiniteWordPath_square_grid_callable(WordDatatype_callable, FiniteWordPath_square_grid, FiniteWord_class):
    pass


# #### Unknown length paths on square grid (experimental) ####

# class WordPath_square_grid_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_square_grid, Word_class):
#    pass

# #### Finite paths on triangle grid ####

class FiniteWordPath_triangle_grid_list(WordDatatype_list, FiniteWordPath_triangle_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='triangle')
        sage: p = P(['a','b','b']); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_triangle_grid_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_triangle_grid_str(WordDatatype_str, FiniteWordPath_triangle_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='triangle')
        sage: p = P('abb'); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_triangle_grid_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_triangle_grid_tuple(WordDatatype_tuple, FiniteWordPath_triangle_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='triangle')
        sage: p = P(('a','b','b')); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_triangle_grid_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_triangle_grid_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_triangle_grid, FiniteWord_class):
    pass


class FiniteWordPath_triangle_grid_iter(WordDatatype_iter, FiniteWordPath_triangle_grid, FiniteWord_class):
    pass


class FiniteWordPath_triangle_grid_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_triangle_grid, FiniteWord_class):
    pass


class FiniteWordPath_triangle_grid_callable(WordDatatype_callable, FiniteWordPath_triangle_grid, FiniteWord_class):
    pass


# #### Finite paths on hexagonal grid ####

class FiniteWordPath_hexagonal_grid_list(WordDatatype_list, FiniteWordPath_hexagonal_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='hexagon')
        sage: p = P(['a','b','b']); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_hexagonal_grid_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_hexagonal_grid_str(WordDatatype_str, FiniteWordPath_hexagonal_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='hexagon')
        sage: p = P('abb'); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_hexagonal_grid_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_hexagonal_grid_tuple(WordDatatype_tuple, FiniteWordPath_hexagonal_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='hexagon')
        sage: p = P(('a','b','b')); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_hexagonal_grid_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_hexagonal_grid_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_hexagonal_grid, FiniteWord_class):
    pass


class FiniteWordPath_hexagonal_grid_iter(WordDatatype_iter, FiniteWordPath_hexagonal_grid, FiniteWord_class):
    pass


class FiniteWordPath_hexagonal_grid_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_hexagonal_grid, FiniteWord_class):
    pass


class FiniteWordPath_hexagonal_grid_callable(WordDatatype_callable, FiniteWordPath_hexagonal_grid, FiniteWord_class):
    pass


# #### Finite paths on cube grid ####

class FiniteWordPath_cube_grid_list(WordDatatype_list, FiniteWordPath_cube_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='cube')
        sage: p = P(['a','b','b']); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_cube_grid_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_cube_grid_str(WordDatatype_str, FiniteWordPath_cube_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='cube')
        sage: p = P('abb'); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_cube_grid_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_cube_grid_tuple(WordDatatype_tuple, FiniteWordPath_cube_grid, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('abcdef', steps='cube')
        sage: p = P(('a','b','b')); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_cube_grid_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_cube_grid_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_cube_grid, FiniteWord_class):
    pass


class FiniteWordPath_cube_grid_iter(WordDatatype_iter, FiniteWordPath_cube_grid, FiniteWord_class):
    pass


class FiniteWordPath_cube_grid_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_cube_grid, FiniteWord_class):
    pass


class FiniteWordPath_cube_grid_callable(WordDatatype_callable, FiniteWordPath_cube_grid, FiniteWord_class):
    pass


# #### Finite paths on north_east ####

class FiniteWordPath_north_east_list(WordDatatype_list, FiniteWordPath_north_east, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('ab', steps='ne')
        sage: p = P(['a','b','b']); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_north_east_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_north_east_str(WordDatatype_str, FiniteWordPath_north_east, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('ab', steps='ne')
        sage: p = P('abb'); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_north_east_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_north_east_tuple(WordDatatype_tuple, FiniteWordPath_north_east, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('ab', steps='ne')
        sage: p = P(('a','b','b')); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_north_east_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_north_east_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_north_east, FiniteWord_class):
    pass


class FiniteWordPath_north_east_iter(WordDatatype_iter, FiniteWordPath_north_east, FiniteWord_class):
    pass


class FiniteWordPath_north_east_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_north_east, FiniteWord_class):
    pass


class FiniteWordPath_north_east_callable(WordDatatype_callable, FiniteWordPath_north_east, FiniteWord_class):
    pass


# #### Finite paths on dyck ####

class FiniteWordPath_dyck_list(WordDatatype_list, FiniteWordPath_dyck, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('ab', steps='dyck')
        sage: p = P(['a','b','b']); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_dyck_list'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_dyck_str(WordDatatype_str, FiniteWordPath_dyck, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('ab', steps='dyck')
        sage: p = P('abb'); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_dyck_str'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_dyck_tuple(WordDatatype_tuple, FiniteWordPath_dyck, FiniteWord_class):
    r"""
    TESTS::

        sage: P = WordPaths('ab', steps='dyck')
        sage: p = P(('a','b','b')); p
        Path: abb
        sage: type(p)
        <class 'sage.combinat.words.paths.FiniteWordPath_dyck_tuple'>
        sage: p == loads(dumps(p))
        True
    """
    pass


class FiniteWordPath_dyck_iter_with_caching(WordDatatype_iter_with_caching, FiniteWordPath_dyck, FiniteWord_class):
    pass


class FiniteWordPath_dyck_iter(WordDatatype_iter, FiniteWordPath_dyck, FiniteWord_class):
    pass


class FiniteWordPath_dyck_callable_with_caching(WordDatatype_callable_with_caching, FiniteWordPath_dyck, FiniteWord_class):
    pass


class FiniteWordPath_dyck_callable(WordDatatype_callable, FiniteWordPath_dyck, FiniteWord_class):
    pass
