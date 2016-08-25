r"""
FindStat - the Combinatorial Statistic Finder.

The FindStat database can be found at http://www.findstat.org .

Fix the following three notions:

- A *combinatorial collection* is a set `S` with interesting combinatorial properties,
- a *combinatorial map* is a combinatorially interesting map `f: S \to S'` between combinatorial collections, and
- a *combinatorial statistic* is a combinatorially interesting map `s: S \to \ZZ`.

You can use the sage interface to FindStat to:

- identify a combinatorial statistic from the values on a few small objects,
- obtain more terms, formulae, references, etc. for a given statistic,
- edit statistics and submit new statistics.

To access the database, use :class:`findstat<FindStat>`::

    sage: findstat
    The Combinatorial Statistic Finder (http://www.findstat.org/)

AUTHORS:

- Martin Rubey (2015): initial version.

A guided tour
-------------

Retrieving information
^^^^^^^^^^^^^^^^^^^^^^

The most straightforward application of the FindStat interface is to
gather information about a combinatorial statistic.  To do this, we
supply :class:`findstat<FindStat>` with a list of `(object, value)`
pairs.  For example::

    sage: PM8 = PerfectMatchings(8)
    sage: r = findstat([(m, m.number_of_nestings()) for m in PM8]); r           # optional -- internet,random
    0: (St000041: The number of nestings of a perfect matching. , [], 105)
    ...

The result of this query is a list (presented as a
:class:`sage.databases.oeis.FancyTuple`) of triples.  The first
element of each triple is a :class:`FindStatStatistic` `s: S \to
\ZZ`, the second element a list of :class:`FindStatMap`'s `f_i: S_i
\to S_{i+1}`, and the third element is an integer::

    sage: (s, list_f, quality) = r[0]                                           # optional -- internet

The precise meaning of the result is as follows:

    The composition `f_n \circ ... \circ f_2 \circ f_1` applied to
    the objects sent to FindStat agrees with `quality` many `(object,
    value)` pairs of `s` in the database.  Moreover, there are no
    other `(object, value)` pairs of `s` stored in the database,
    i.e., there is no disagreement of values.

Put differently, if `quality` is not too small it is likely that the
statistic sent to FindStat equals `s \circ f_n \circ ... \circ f_2 \circ f_1`.

In the case at hand, the list of maps is empty and the integer
`quality` equals the number of `(object, value)` pairs passed to
FindStat.  This means, that the set of `(object, value)` pairs of the
statistic `s` as stored in the FindStat database is a superset of the
data sent.  We can now retrieve the description from the database::

    sage: print(s.description())                                         # optional -- internet,random
    The number of nestings of a perfect matching.
    <BLANKLINE>
    <BLANKLINE>
    This is the number of pairs of edges $((a,b), (c,d))$ such that $a\le c\le d\le b$. i.e., the edge $(c,d)$ is nested inside $(a,b)$.

and check the references::

    sage: s.references()                                                        # optional -- internet,random
    0: [1] [[MathSciNet:1288802]]
    1: [2] [[MathSciNet:1418763]]

If you prefer, you can look at this information also in your browser::

    sage: findstat(41).browse()                                                 # optional -- webbrowser

Another interesting possibility is to look for equidistributed
statistics.  Instead of submitting a list of pairs, we pass a pair of
lists::

    sage: r = findstat((PM8, [m.number_of_nestings() for m in PM8])); r         # optional -- internet,random
    0: (St000041: The number of nestings of a perfect matching. , [], 105)
    1: (St000042: The number of crossings of a perfect matching. , [], 105)
    ...

This results tells us that the database contains another entriy that is
equidistributed with the number of nestings on perfect matchings of
length `8`, namely the number of crossings.

Let us now look at a slightly more complicated example, where the
submitted statistic is the composition of a sequence of combinatorial
maps and a statistic known to FindStat.  We use the occasion to
advertise yet another way to pass values to FindStat::

    sage: r = findstat(Permutations(4), lambda pi: pi.saliances()[0]); r        # optional -- internet,random
    0: (St000051: The size of the left subtree. , [Mp00069: complement, Mp00061: to increasing tree], 24)
    ...
    sage: (s, list_f, quality) = r[0]                                           # optional -- internet

To obtain the value of the statistic sent to FindStat on a given
object, apply the maps in the list in the given order to this object,
and evaluate the statistic on the result.  For example, let us check
that the result given by FindStat agrees with our statistic on the
following permutation::

    sage: pi = Permutation([3,1,4,5,2]); pi.saliances()[0]
    3

We first have to find out, what the maps and the statistic actually do::

    sage: print(s.description())                                            # optional -- internet,random
    The size of the left subtree.

    sage: print(s.code())                                                 # optional -- internet,random
    def statistic(T):
        return T[0].node_number()

    sage: print(list_f[0].code() + "\r\n" + list_f[1].code())               # optional -- internet,random
    def complement(elt):
        n = len(elt)
        return elt.__class__(elt.parent(), map(lambda x: n - x + 1, elt) )
    <BLANKLINE>
    def increasing_tree_shape(elt, compare=min):
        return elt.increasing_tree(compare).shape()

So, the following should coincide with what we sent FindStat::

    sage: pi.complement().increasing_tree_shape()[0].node_number()
    3

Editing and submitting statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Of course, often a statistic will not be in the database::

    sage: findstat([(d, randint(1,1000)) for d in DyckWords(4)])                # optional -- internet
    a new statistic on Cc0005: Dyck paths

In this case, and if the statistic might be "interesting", please
consider submitting it to the database using
:meth:`FindStatStatistic.submit`.

Also, you may notice omissions, typos or even mistakes in the
description, the code and the references.  In this case, simply
replace the value by using :meth:`FindStatStatistic.set_description`,
:meth:`FindStatStatistic.set_code` or
:meth:`FindStatStatistic.set_references`, and then
:meth:`FindStatStatistic.submit` your changes for review by the
FindStat team.

Classes and methods
-------------------
"""
#*****************************************************************************
#       Copyright (C) 2015 Martin Rubey <martin.rubey@tuwien.ac.at>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute

from sage.categories.sets_cat import Sets
from sage.structure.sage_object import SageObject

from sage.misc.misc import verbose
from sage.rings.integer import Integer
from sage.databases.oeis import FancyTuple

from ast import literal_eval
from collections import OrderedDict
import re
import webbrowser
import tempfile
import time
import inspect
import json
import cgi

# import compatible with py2 and py3
from six.moves.urllib.parse import urlencode
from six.moves.urllib.request import Request, urlopen
from six.moves.urllib.error import HTTPError

# Combinatorial collections
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix, AlternatingSignMatrices
from sage.combinat.binary_tree import BinaryTree, BinaryTrees
from sage.combinat.core import Core, Cores
from sage.combinat.dyck_word import DyckWord, DyckWords
from sage.combinat.root_system.cartan_type import CartanType_abstract, CartanType
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern, GelfandTsetlinPatterns
from sage.graphs.graph import Graph
from sage.combinat.composition import Composition, Compositions
from sage.combinat.partition import Partition, Partitions
from sage.combinat.ordered_tree import OrderedTree, OrderedTrees
from sage.combinat.parking_functions import ParkingFunction, ParkingFunction_class, ParkingFunctions
from sage.combinat.perfect_matching import PerfectMatching, PerfectMatchings
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.posets.posets import Poset, FinitePoset
from sage.combinat.posets.poset_examples import posets
from sage.combinat.tableau import SemistandardTableau, SemistandardTableaux, StandardTableau, StandardTableaux
from sage.combinat.set_partition import SetPartition, SetPartitions
from sage.graphs.graph_generators import graphs
from sage.combinat.words.word import Word
from sage.combinat.words.words import Words
from sage.combinat.words.abstract_word import Word_class

######################################################################
# the FindStat URLs
FINDSTAT_URL                        = 'http://www.findstat.org/'
FINDSTAT_URL_RESULT                 = FINDSTAT_URL + "StatisticFinder/Result/"
FINDSTAT_URL_LOGIN                  = FINDSTAT_URL + "StatisticFinder?action=login"
FINDSTAT_URL_NEW                    = FINDSTAT_URL + 'StatisticsDatabase/NewStatistic/'
FINDSTAT_URL_EDIT                   = FINDSTAT_URL + 'StatisticsDatabase/EditStatistic/'
FINDSTAT_URL_BROWSE                 = FINDSTAT_URL + 'StatisticsDatabase/'

FINDSTAT_URL_DOWNLOADS              = 'http://downloads.findstat.org/'
FINDSTAT_URL_DOWNLOADS_STATISTICS   = FINDSTAT_URL_DOWNLOADS + "statistics/%s.json"
FINDSTAT_URL_DOWNLOADS_COLLECTIONS  = FINDSTAT_URL_DOWNLOADS + "collections.json"
FINDSTAT_URL_DOWNLOADS_MAPS         = FINDSTAT_URL_DOWNLOADS + "maps.json"

######################################################################
# the number of values FindStat allows to search for at most
FINDSTAT_MAX_VALUES = 200
# the number of values FindStat needs at least to search for
FINDSTAT_MIN_VALUES = 3
# the number of maps that FindStat should compose at most to find a match
FINDSTAT_MAX_DEPTH = 5
# the number of values FindStat allows to submit at most
FINDSTAT_MAX_SUBMISSION_VALUES = 1200

# the fields of the FindStat database we expect
FINDSTAT_STATISTIC_IDENTIFIER                   = 'StatisticIdentifier'
FINDSTAT_STATISTIC_COLLECTION                   = 'StatisticCollection'
FINDSTAT_STATISTIC_DATA                         = 'StatisticData'
FINDSTAT_STATISTIC_GENERATING_FUNCTION          = 'StatisticGeneratingFunction'
FINDSTAT_STATISTIC_NAME                        = 'StatisticTitle'
FINDSTAT_STATISTIC_DESCRIPTION                  = 'StatisticDescription'
FINDSTAT_STATISTIC_REFERENCES                   = 'StatisticReferences'
FINDSTAT_STATISTIC_CODE                         = 'StatisticCode'
FINDSTAT_STATISTIC_ORIGINAL_AUTHOR              = 'StatisticOriginalAuthor' # unused, designates a dictionary with Name, Time
FINDSTAT_STATISTIC_UPDATE_AUTHOR                = 'StatisticUpdateAuthor'   # unused, designates a dictionary with Name, Time

FINDSTAT_POST_AUTHOR                            = 'StatisticAuthor' # designates the name of the author
FINDSTAT_POST_EMAIL                             = 'StatisticEmail'
FINDSTAT_POST_SAGE_CELL                         = 'SageCellField'   # currently only used as post key
FINDSTAT_POST_EDIT                              = 'EDIT'            # only used as post key

FINDSTAT_COLLECTION_IDENTIFIER                  = 'CollectionIdentifier'
FINDSTAT_COLLECTION_NAME                        = 'CollectionName'
FINDSTAT_COLLECTION_NAME_PLURAL                 = 'CollectionNamePlural'
FINDSTAT_COLLECTION_NAME_WIKI                   = 'CollectionNameWiki'
FINDSTAT_COLLECTION_LEVELS                      = 'CollectionLevelsSizes'

FINDSTAT_MAP_IDENTIFIER                         = 'MapIdentifier'
FINDSTAT_MAP_NAME                               = 'MapName'
FINDSTAT_MAP_DESCRIPTION                        = 'MapDescription'
FINDSTAT_MAP_DOMAIN                             = 'MapDomain'
FINDSTAT_MAP_CODOMAIN                           = 'MapCodomain'
FINDSTAT_MAP_CODE                               = 'MapCode'
FINDSTAT_MAP_CODE_NAME                          = 'MapSageName'

FINDSTAT_QUERY_MATCHES                          = 'QueryMatches'
FINDSTAT_QUERY_MATCHING_DATA                    = 'QueryMatchingData'
FINDSTAT_QUERY_MAPS                             = 'QueryMaps'

# the entries of this list are required as post arguments for submitting or editing a statistic
FINDSTAT_EDIT_FIELDS = set([FINDSTAT_STATISTIC_IDENTIFIER,
                            FINDSTAT_STATISTIC_COLLECTION,
                            FINDSTAT_STATISTIC_DATA,
                            FINDSTAT_STATISTIC_DESCRIPTION,
                            FINDSTAT_STATISTIC_REFERENCES,
                            FINDSTAT_STATISTIC_CODE,
                            FINDSTAT_POST_AUTHOR,
                            FINDSTAT_POST_EMAIL,
                            FINDSTAT_POST_SAGE_CELL,
                            FINDSTAT_POST_EDIT])

# separates name from description
FINDSTAT_SEPARATOR_NAME = "\r\n"
# separates references
FINDSTAT_SEPARATOR_REFERENCES = "\r\n"

######################################################################

# the format string for using POST
# WARNING: we use cgi.escape to avoid injection problems, thus we expect double quotes as field delimiters.
FINDSTAT_POST_HEADER = """
<script src="http://www.google.com/jsapi"></script>
<script>
    google.load("jquery", "1.3.2");
</script>

<script>
    $(document).ready(function() {$("#form").submit(); });
</script>
"""
FINDSTAT_NEWSTATISTIC_FORM_HEADER = '<form id="form" name="NewStatistic" action="%s" enctype="multipart/form-data" method="post" />'
FINDSTAT_NEWSTATISTIC_FORM_FORMAT = '<input type="hidden" name="%s" value="%s" />'
FINDSTAT_NEWSTATISTIC_FORM_FOOTER = '</form>'

######################################################################
class FindStat(SageObject):
    r"""
    The Combinatorial Statistic Finder.

    :class:`FindStat` is a class representing results of queries to
    the FindStat database.  This class is also the entry point to
    edit statistics and new submissions.  Use the shorthand
    :class:`findstat<FindStat>` to call it.

    INPUT:

    One of the following:

    - an integer or a string representing a valid FindStat identifier
      (e.g. 45 or 'St000045').  The keyword arguments ``depth`` and
      ``max_values`` are ignored.

    - a list of pairs of the form (object, value), or a dictionary
      from sage objects to integer values.  The keyword arguments
      ``depth`` and ``max_values`` are passed to the finder.

    - a list of pairs of the form (list of objects, list of values),
      or a single pair of the form (list of objects, list of values).
      In each pair there should be as many objects as values.  The
      keyword arguments ``depth`` and ``max_values`` are passed to
      the finder.

    - a collection and a list of pairs of the form (string, value),
      or a dictionary from strings to integer values.  The keyword
      arguments ``depth`` and ``max_values`` are passed to the
      finder.  This should only be used if the collection is not yet
      supported.

    - a collection and a list of pairs of the form (list of strings,
      list of values), or a single pair of the form (list of strings,
      list of values).  In each pair there should be as many strings
      as values.  The keyword arguments ``depth`` and ``max_values``
      are passed to the finder.  This should only be used if the
      collection is not yet supported.

    - a collection and a callable.  The callable is used to generate
      ``max_values`` (object, value) pairs.  The number of terms
      generated may also be controlled by passing an iterable
      collection, such as ``Permutations(3)``.  The keyword arguments
      ``depth`` and ``max_values`` are passed to the finder.

    OUTPUT:

    An instance of a :class:`FindStatStatistic`, represented by

    - the FindStat identifier together with its name, or

    - a list of triples, each consisting of

        - the statistic

        - a list of strings naming certain maps

        - a number which says how many of the values submitted agree
          with the values in the database, when applying the maps in
          the given order to the object and then computing the
          statistic on the result.

    EXAMPLES:

    A particular statistic can be retrieved by its St-identifier or
    number::

        sage: findstat('St000041')                                              # optional -- internet,random
        St000041: The number of nestings of a perfect matching.

        sage: findstat(51)                                                      # optional -- internet,random
        St000051: The size of the left subtree.

    The database can be searched by providing a list of pairs::

        sage: q = findstat([(pi, pi.length()) for pi in Permutations(4)]); q    # optional -- internet,random
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
        1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)
        ...

    or a dictionary::

        sage: p = findstat({pi: pi.length() for pi in Permutations(4)}); p      # optional -- internet,random
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
        1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)
        ...

    Note however, that the results of these two queries are not
    necessarily the same, because we compare queries by the data
    sent, and the ordering of the data might be different::

        sage: p == q                                                            # optional -- internet
        False

    Another possibility is to send a collection and a function.  In
    this case, the function is applied to the first few objects of
    the collection::

        sage: findstat("Permutations", lambda pi: pi.length())                  # optional -- internet,random
        0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 200)
        ...

    To search for a distribution, send a list of lists, or a single pair::

        sage: S4 = Permutations(4); findstat((S4, [pi.length() for pi in S4]))  # optional -- internet,random
        0: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [], 24)
        1: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
        ...

    Note that there is a limit, ``FINDSTAT_MAX_DEPTH``, on the number
    of elements that may be submitted to FindStat, which is currently
    200.  Therefore, the interface tries to truncate queries
    appropriately, but this may be impossible, especially with
    distribution searches::

        sage: S6 = Permutations(6); S6.cardinality()                            # optional -- internet
        720
        sage: findstat((S6, [1 for a in S6]))                                   # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: After discarding elements not in the range, too few (=0) values remained to send to FindStat.
    """
    def __init__(self):
        r"""
        Initialize the database.

        In particular, we set up a cache from integers to
        :class:`FindStatStatistic` instances that avoids retrieving
        the same statistic over and over again.

        TESTS::

            sage: findstat
            The Combinatorial Statistic Finder (http://www.findstat.org/)
        """
        self._statistic_cache = dict()

        # user credentials if provided
        self._user_name  = ""
        self._user_email = ""

    def __call__(self, query_1, query_2=None, depth=2, max_values=FINDSTAT_MAX_VALUES):
        r"""
        Return an instance of a :class:`FindStatStatistic`.

        This should be the only way to access
        :class:`FindStatStatistic`.  We do the preprocessing and
        (some) checking of the data here.  Then we call the
        appropriate method of :class:`FindStatStatistic` to launch
        the query.  Thus, in principle :class:`FindStatStatistic`
        could be called if no checking is desired.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: findstat(FindStatCollection("Permutations"), lambda pi: pi.length())    # optional -- internet,random
            0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 200)
            1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 200)
            2: (St000067: The inversion number of the alternating sign matrix., [Mp00063: to alternating sign matrix], 152)
            3: (St000246: The [[/Permutations/Inversions|number of non-inversions]] of a permutation., [Mp00064: reverse], 200)
            4: (St000008: The major index of the composition., [Mp00062: inversion-number to major-index bijection, Mp00071: descent composition], 200)

        TESTS::

            sage: findstat("Permutations", lambda x: 1, depth="x")
            Traceback (most recent call last):
            ...
            ValueError: The depth of a FindStat query must be a non-negative integer less than or equal to 5.

            sage: findstat("Permutations", lambda x: 1, depth=100)
            Traceback (most recent call last):
            ...
            ValueError: The depth of a FindStat query must be a non-negative integer less than or equal to 5.

            sage: findstat("Permutations", 1)                                   # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: The given arguments, Permutations and 1, cannot be used for a FindStat search.

            sage: S = Permutation
            sage: findstat([[S([1,2,3]),S([1,3,2])],[1,1]])                     # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: After discarding elements not in the range, too few (=2) values remained to send to FindStat.

            sage: findstat(([S([1,3,2]), S([1,2,3]), S([1,3,2])], [1,1,1]))     # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: FindStat expects that every object occurs at most once.

            sage: findstat([(S([1,2]), 1), ([S([1,3,2]), S([1,2])], [2,3])])    # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: FindStat expects that every object occurs at most once.
        """
        try:
            depth = int(depth)
            assert 0 <= depth <= FINDSTAT_MAX_DEPTH
        except (ValueError, AssertionError):
            raise ValueError("The depth of a FindStat query must be a non-negative integer less than or equal to %i." %FINDSTAT_MAX_DEPTH)

        def get_collection(collection=None, element=None):
            if collection is None:
                collection = FindStatCollection(element)
                return (collection, collection.to_string())
            else:
                return (FindStatCollection(collection), lambda x: x)

        def query_by_dict(query, collection=None):
            # we expect a dictionary from objects or strings to
            # integers
            l = query.iteritems()
            (key, value) = next(l)

            (collection, to_str) = get_collection(collection, key)

            data = ([([key], [to_str(key)], [Integer(value)])] +
                    [([key], [to_str(key)], [Integer(value)]) for (key, value) in l])

            first_terms = [(key[0], value[0]) for (key, key_str, value) in data]

            return FindStatStatistic(id=0, data=data,
                                     first_terms=first_terms,
                                     collection=collection,
                                     depth=depth)._find_by_values(max_values=max_values)

        def query_by_iterable(query, collection=None):
            # either a pair (list of objects, list of integers)
            # or a list of such or (object, integer) pairs

            # values must always be converted to lists because
            # otherwise we get a trailing comma when printing
            if (len(query) == 2 and
                isinstance(query[1], (list, tuple)) and
                len(query[1]) != 0 and
                isinstance(query[1][0], (int, Integer))):
                # just a single pair, i.e., a pure distribution query

                (collection, to_str) = get_collection(collection, query[0][0])
                data = [(query[0], map(to_str, query[0]), map(Integer, query[1]))]
                if len(data[0][0]) != len(data[0][2]):
                    raise ValueError("FindStat expects the same number of objects as values.")
                if len(set(data[0][1])) != len(data[0][2]):
                    raise ValueError("FindStat expects that every object occurs at most once.")

                return FindStatStatistic(id=0, data=data,
                                         collection=collection,
                                         depth=depth)._find_by_values(max_values=max_values)
            else:
                (key, value) = query[0]
                if isinstance(value, (list, tuple)):
                    (collection, to_str) = get_collection(collection, key[0])
                else:
                    (collection, to_str) = get_collection(collection, key)

                data = []
                is_statistic = True
                for (key, value) in query:
                    if isinstance(value, (list, tuple)):
                        if len(key) != len(value):
                            raise ValueError("FindStat expects the same number of objects as values.")
                        if len(value) != 1:
                            is_statistic = False
                        data += [(key, map(to_str, key), map(Integer, value))]
                    else:
                        data += [([key], [to_str(key)], [Integer(value)])]

                all_elements = [e for (elements, elements_str, values) in data for e in elements_str]
                if len(set(all_elements)) != len(all_elements):
                    raise ValueError("FindStat expects that every object occurs at most once.")

                if is_statistic:
                    return FindStatStatistic(id=0, data=data,
                                             collection=collection,
                                             first_terms=query,
                                             depth=depth)._find_by_values(max_values=max_values)
                else:
                    return FindStatStatistic(id=0, data=data,
                                             collection=collection,
                                             depth=depth)._find_by_values(max_values=max_values)

        if query_2 is None:
            if isinstance(query_1, str):
                if re.match('^St[0-9]{6}$', query_1):
                    return self._statistic_find_by_id_cached(Integer(query_1[2:].lstrip("0")))
                else:
                    raise ValueError("The value %s is not a valid FindStat statistic identifier." %query_1)

            elif isinstance(query_1, (int, Integer)):
                return self._statistic_find_by_id_cached(query_1)

            elif isinstance(query_1, dict):
                return query_by_dict(query_1)

            elif isinstance(query_1, (list, tuple)):
                return query_by_iterable(query_1)

            else:
                raise ValueError("The given query, %s, cannot be used for a FindStat search." %query_1)

        else:
            (collection, to_str) = get_collection(None, query_1)
            if isinstance(query_2, dict):
                return query_by_dict(query_2, collection)

            elif isinstance(query_2, (list, tuple)):
                return query_by_iterable(query_2, collection)

            elif callable(query_2):
                first_terms = collection.first_terms(query_2, max_values=max_values)
                data = [([key], [to_str(key)], [Integer(value)]) for (key, value) in first_terms]
                try:
                    code = inspect.getsource(query_2)
                except IOError:
                    _ = verbose("inspect.getsource could not get code from function provided", caller_name='FindStat')
                    code = ""
                return FindStatStatistic(id=0, first_terms=first_terms,
                                         data=data, function=query_2, code=code,
                                         collection=collection,
                                         depth=depth)._find_by_values(max_values=max_values)
            else:
                raise ValueError("The given arguments, %s and %s, cannot be used for a FindStat search." %(query_1, query_2))

    def __repr__(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: findstat
            The Combinatorial Statistic Finder (http://www.findstat.org/)
        """
        return "The Combinatorial Statistic Finder (%s)" % FINDSTAT_URL

    def browse(self):
        r"""
        Open the FindStat web page in a browser.

        EXAMPLES::

            sage: findstat.browse()                                             # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL)

    def set_user(self, name=None, email=None):
        r"""
        Set the user for this session.

        INPUT:

        - ``name`` -- the name of the user.

        - ``email`` -- an email address of the user.

        This information is used when submitting a statistic with
        :meth:`FindStatStatistic.submit`.

        EXAMPLES::

            sage: findstat.set_user(name="Anonymous", email="invalid@org")

        .. NOTE::

            It is usually more convenient to login into the FindStat
            web page using the :meth:`login` method.
        """
        if not isinstance(name, str):
            raise ValueError("The given name is not a string.")
        if not isinstance(email, str):
            raise ValueError("The given email address is not a string.")
        self._user_name  = name
        self._user_email = email

    def login(self):
        r"""
        Open the FindStat login page in a browser.

        EXAMPLES::

            sage: findstat.login()                                              # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_LOGIN)

    ######################################################################

    def _statistic_find_by_id_cached(self, id):
        r"""
        INPUT:

        - ``id`` -- an integer designating the FindStat id of a statistic.

        OUTPUT:

        An instance of :class:`FindStatStatistic`.

        .. TODO::

            this method caches the statistics.  It may make sense to
            provide a method that clears the cache, or reloads a
            single statistic.

        TESTS::

            sage: findstat(41).set_description("")                              # optional -- internet, indirect doctest
            sage: findstat(41).description() == ""                              # optional -- internet, indirect doctest
            True
        """
        if id > 0:
            if id not in self._statistic_cache.keys():
                self._statistic_cache[id] = FindStatStatistic(id)._find_by_id()
            return self._statistic_cache[id]
        else:
            raise ValueError("A FindStat statistic identifier must be at least 1.")

######################################################################

class FindStatStatistic(SageObject):
    r"""
    The class of FindStat statistics.

    Do not instantiate this class directly.  Instead, use
    :class:`findstat<FindStat>`.
    """
    def __init__(self, id, first_terms=None, data=None, function=None, code="", collection=None, depth=None):
        r"""
        Initialize a FindStat query for a statistic from preprocessed
        data.

        INPUT:

        - ``id`` -- an integer designating the FindStat id of a
          statistic, or 0.

        - ``first_terms`` -- (optional) a list of (object, value)
          pairs, see :meth:`first_terms`.

        - ``data`` -- (optional) a list of pairs of the form (list of
          objects represented as strings, list of values), see
          :meth:`data`.

        - ``function`` -- (optional) a function taking a sage object
          as input and returning the value of the statistic on this
          object, see :meth:`function`.

        - ``code`` -- (optional) a string containing code (possibly
          pseudocode or code for a different computer algebra
          system), see :meth:`code`.

        - ``collection`` -- (optional) an instance of
          :class:`FindStatCollection`, see :meth:`collection`.

        - ``depth`` -- (optional) an integer between 0 and
          FINDSTAT_MAX_DEPTH, which determines how many maps FindStat
          should compose at most to find a match.

        This method by itself does not launch the query, because
        there are two rather different queries possible, see
        :meth:`_find_by_id` and :meth:`_find_by_values`.

        TESTS::

            sage: from sage.databases.findstat import FindStatStatistic
            sage: FindStatStatistic(1)._find_by_id()                            # optional -- internet
            St000001: The number of ways to write a permutation as a minimal length product of simple transpositions.
        """
        self._depth = depth
        self._query = None
        self._modified = False

        self._id = id
        self._result = None

        self._first_terms = first_terms
        self._data = data
        self._function = function
        self._code = code
        self._collection = collection

        self._description = ""
        self._references = ""

    def __repr__(self):
        r"""
        Return the representation of the FindStat query.

        OUTPUT:

        A string.  If the query was by identifier, the identifier and
        the name of the statistic.  If the statistic was modified
        (see :meth:`modified`) this is also indicated.

        If the query was by values and at least one match was found,
        the list of matches.

        Otherwise, if no match was found, a message saying this.

        EXAMPLES::

            sage: findstat(1)                                                   # optional -- internet
            St000001: ...

            sage: findstat([(pi, randint(1,100)) for pi in Permutations(3)])    # optional -- internet
            a new statistic on Cc0001: Permutations

            sage: findstat([(pi, pi(1)) for pi in Permutations(3)])             # optional -- internet
            0: (St000054: ...
            1: ...
            2: ...

        """
        if self._query == "ID":
            if self._modified:
                return "%s(modified): %s" % (self.id_str(), self.name())
            else:
                return "%s: %s" % (self.id_str(), self.name())

        elif self._query == "data":
            if len(self._result) == 0:
                return "a new statistic on " + repr(self._collection)
            else:
                return repr(self._result)

        else:
            raise ValueError("FindStatStatistic._query should be either 'ID' or 'data', but is %s.  This should not happen.  Please send an email to the developers." %self._query)

    def __eq__(self, other):
        """Return ``True`` if ``self`` is equal to ``other`` and ``False``
        otherwise.

        INPUT:

        - ``other`` -- a FindStat query, i.e., instance of :class:`FindStatStatistic`.

        OUTPUT:

        A boolean.

        Two queries are considered equal if all of the following
        applies:

        - the queries are both of the same type, i.e., both are by
          identifier or both are by values,

        - if the queries are by identifier, the identifiers agree and
          the statistics are both unmodified,

        - if the queries are by values, the submitted data are the
          same and the statistic data are both unmodified.

        .. TODO::

            this is *very* rudimentary

        EXAMPLES::

            sage: findstat(1) == findstat(41)                                   # optional -- internet
            False

            sage: r1 = findstat(Permutations(3), lambda pi: pi.saliances()[0])  # optional -- internet
            sage: r2 = findstat(Permutations(4), lambda pi: pi.saliances()[0])  # optional -- internet
            sage: r1 == r2                                                      # optional -- internet
            False
        """
        if (not isinstance(other, FindStatStatistic)):
            return False
        if self._query == "ID" and other._query == "ID":
            if self._modified or other._modified:
                return False
            else:
                return self._id == other._id
        elif self._query == "data" and other._query == "data":
            if self._modified or other._modified:
                return False
            else:
                return self._data == other._data
        else:
            return False

    def __ne__(self, other):
        """Determine whether ``other`` is a different query.

        INPUT:

        - ``other`` -- a FindStat query, i.e., instance of
          :class:`FindStatStatistic`..

        OUTPUT:

        A boolean.

        SEEALSO:

        :meth:`__eq__`

        EXAMPLES::

            sage: r1 = findstat(Permutations(3), lambda pi: pi.saliances()[0])  # optional -- internet
            sage: r2 = findstat(Permutations(4), lambda pi: pi.saliances()[0])  # optional -- internet
            sage: r1 != r2                                                      # optional -- internet
            True

        """
        return not self == other


    ######################################################################
    # query = "ID"
    ######################################################################

    def _find_by_id(self):
        r"""
        Retrieve the statistic matching ``self._id``.

        OUTPUT:

        The statistic ``self``.

        Expects that ``_id`` is a valid identifier.  Overwrites all
        variables associated with the statistic, such as
        ``_description``, ``_code``, ``_references``, etc.

        TESTS::

            sage: findstat(999999)                                              # optional -- internet, indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: St999999 is not a FindStat statistic identifier.
        """
        self._query = "ID"

        # get the database entry from FindStat
        url = FINDSTAT_URL_DOWNLOADS_STATISTICS %self.id_str()
        _ = verbose("Fetching URL %s ..." %url, caller_name='FindStat')
        try:
            self._raw = json.load(urlopen(url), object_pairs_hook=OrderedDict)
        except HTTPError as error:
            if error.code == 404:
                raise ValueError("%s is not a FindStat statistic identifier." %self.id_str())
            else:
                raise

        self._description           = self._raw[FINDSTAT_STATISTIC_DESCRIPTION].encode("utf-8")
        self._name                  = self._raw[FINDSTAT_STATISTIC_NAME].encode("utf-8")
        self._references            = self._raw[FINDSTAT_STATISTIC_REFERENCES].encode("utf-8")
        self._collection            = FindStatCollection(self._raw[FINDSTAT_STATISTIC_COLLECTION])
        self._code                  = self._raw[FINDSTAT_STATISTIC_CODE]

        from ast import literal_eval
        gf                          = self._raw[FINDSTAT_STATISTIC_GENERATING_FUNCTION]
        self._generating_functions_dict  = { literal_eval(key):
                                             { literal_eval(inner_key): inner_value
                                               for inner_key, inner_value in value.iteritems() }
                                             for key, value in gf.iteritems() }

        from_str = self._collection.from_string()
        # we want to keep FindStat's ordering here!
        self._first_terms = [(from_str(obj), Integer(val)) for (obj, val) in self._raw[FINDSTAT_STATISTIC_DATA].iteritems()]
        return self

    ######################################################################
    # query = "data"
    ######################################################################

    def _find_by_values(self, max_values=FINDSTAT_MAX_VALUES):
        r"""
        Retrieve the statistics matching ``self._data``.

        OUTPUT:

        The query ``self``.

        Expects that ``_data`` is a list of triples of the form (list
        of objects, list of objects represented as strings, list of
        values), each containing as many values as objects, and that
        ``_collection`` is appropriately set.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: from sage.databases.findstat import FindStatStatistic
            sage: collection = FindStatCollection("Dyck paths")                                     # optional -- internet
            sage: to_str = collection.to_string()                                                   # optional -- internet
            sage: query = {dw:dw.area() for dw in DyckWords(4)}
            sage: data = [([key], [to_str(key)], [value]) for (key, value) in query.iteritems()]    # optional -- internet
            sage: first_terms = [(key, value) for (key, value) in query.iteritems()]

            sage: FindStatStatistic(id=0,data=data, first_terms = first_terms, collection = collection, depth=0)._find_by_values() # optional -- internet
            0: (St000012: The area of a Dyck path., [], 14)
        """
        self._query = "data"

        # FindStat allows to search for at most FINDSTAT_MAX_VALUES
        # values.  For the user's convenience, from _data, we take the
        # first min(max_values, FINDSTAT_MAX_VALUES) such that all
        # elements are in the precomputed range
        # the internal _data is not modified
        data = []
        total = min(max_values, FINDSTAT_MAX_VALUES)
        in_range = self._collection.in_range
        in_range_counter = 0
        for (elements, elements_str, values) in self._data:
            if total >= len(elements):
                if all(in_range(e) for e in elements):
                    data                += [(elements, elements_str, values)]
                    in_range_counter    += len(elements)
                    total               -= len(elements)

        if in_range_counter < FINDSTAT_MIN_VALUES:
            raise ValueError("After discarding elements not in the range, too few (=%s) values remained to send to FindStat." %in_range_counter)

        url = FINDSTAT_URL_RESULT + self._collection._url_name + "/"

        stat = [(elements_str, str(values)[1:-1]) for (elements, elements_str, values) in data]

        stat_str = "\n".join(["\n".join(keys) + "\n====> " + values for (keys, values) in stat])
        _ = verbose("Sending the following data to FindStat\r\n %s" %stat_str, caller_name='FindStat')

        values = urlencode({"freedata": stat_str, "depth": str(self._depth), "caller": "Sage"})
        _ = verbose("Fetching URL %s with encoded data %s" %(url, values), caller_name='FindStat')

        request = Request(url, data=values)
        _ = verbose("Requesting %s" %request, caller_name='FindStat')

        response = urlopen(request)
        _ = verbose("Response was %s" %response.info(), caller_name='FindStat')

        try:
            result = json.load(response)
        except Exception as e:
            raise IOError("FindStat did not answer with a json response: %s" %e)

        self._result = FancyTuple((findstat(match[FINDSTAT_STATISTIC_IDENTIFIER]),
                                   [FindStatMap(mp[FINDSTAT_MAP_IDENTIFIER]) for mp in match[FINDSTAT_QUERY_MAPS]],
                                   len(match[FINDSTAT_QUERY_MATCHING_DATA]))
                                  for match in result[FINDSTAT_QUERY_MATCHES])

        return self

    def _raise_error_modifying_statistic_with_perfect_match(self):
        r"""
        Raise an error when there is a result with depth 0.

        TESTS::

            sage: s = findstat([(pi, pi[0]) for pi in Permutations(3)])         # optional -- internet
            sage: s._raise_error_modifying_statistic_with_perfect_match()       # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: Your input data matches St000054.  Consider modifying this statistic instead.

            sage: s = findstat(1)                                               # optional -- internet
            sage: s._raise_error_modifying_statistic_with_perfect_match()       # optional -- internet

            sage: s = findstat([(d, randint(1,1000)) for d in DyckWords(4)])    # optional -- internet
            sage: s._raise_error_modifying_statistic_with_perfect_match()       # optional -- internet
        """
        if self._query == "data" and len(self._result) > 0 and len(self._result[0][1]) == 0:
            raise ValueError("Your input data matches %s.  Consider modifying this statistic instead." %self._result[0][0].id_str())

    ######################################################################

    def __getitem__(self, key):
        r"""
        Return the match corresponding to ``key``.

        INPUT:

        - ``key`` -- an integer.

        OUTPUT:

        The match corresponding to ``key``.  Raise an error if the
        query was by identifier.

        EXAMPLES::

            sage: q = findstat([(pi, pi.length()) for pi in Permutations(4)]); q          # optional -- internet,random
                0: (St000018: The [[/Permutations/Inversions|number of inversions]] of a permutation., [], 24)
                1: (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)
                ...

            sage: q[1]                                                          # optional -- internet,random
            (St000004: The [[/Permutations/Descents-Major|major index]] of a permutation., [Mp00062: inversion-number to major-index bijection], 24)

        """
        if self._query == "ID":
            raise TypeError("Use 'first_terms' to access the values of a FindStatStatistic.")

        elif self._query == "data":
            return self._result[key]

        else:
            raise ValueError("FindStatStatistic._query should be either 'ID' or 'data', but is %s.  This should not happen.  Please send an email to the developers." %self._query)

    def id(self):
        r"""
        Return the FindStat identifier of the statistic.

        OUTPUT:

        The FindStat identifier of the statistic (or 0), as an integer.

        EXAMPLES::

            sage: findstat(1).id()                                              # optional -- internet
            1
        """
        return self._id

    def id_str(self):
        r"""
        Return the FindStat identifier of the statistic.

        OUTPUT:

        The FindStat identifier of the statistic (or 'St000000'), as a string.

        EXAMPLES::

            sage: findstat(1).id_str()                                          # optional -- internet
            'St000001'
        """
        id = str(self._id)
        return 'St000000'[:-len(id)] + id

    def data(self):
        r"""
        Return the data used for querying the FindStat database.

        OUTPUT:

        The data provided by the user to query the FindStat database.
        When the database was searched using an identifier, ``data``
        is ``None``.

        EXAMPLES::

            sage: S4 = Permutations(4); findstat((S4, [pi.length() for pi in S4])).data() # optional -- internet
            [(Standard permutations of 4,
              ['[1, 2, 3, 4]',
               '[1, 2, 4, 3]',
               '[1, 3, 2, 4]',
               '[1, 3, 4, 2]',
               '[1, 4, 2, 3]',
               '[1, 4, 3, 2]',
               '[2, 1, 3, 4]',
               '[2, 1, 4, 3]',
               '[2, 3, 1, 4]',
               '[2, 3, 4, 1]',
               '[2, 4, 1, 3]',
               '[2, 4, 3, 1]',
               '[3, 1, 2, 4]',
               '[3, 1, 4, 2]',
               '[3, 2, 1, 4]',
               '[3, 2, 4, 1]',
               '[3, 4, 1, 2]',
               '[3, 4, 2, 1]',
               '[4, 1, 2, 3]',
               '[4, 1, 3, 2]',
               '[4, 2, 1, 3]',
               '[4, 2, 3, 1]',
               '[4, 3, 1, 2]',
               '[4, 3, 2, 1]'],
              [0, 1, 1, 2, 2, 3, 1, 2, 2, 3, 3, 4, 2, 3, 3, 4, 4, 5, 3, 4, 4, 5, 5, 6])]
        """
        return self._data

    def modified(self):
        r"""
        Return whether the statistic was modified.

        OUTPUT:

        True, if the statistic was modified using
        :meth:`set_description`, :meth:`set_code`,
        :meth:`set_references`, etc.  False otherwise.

        EXAMPLES::

            sage: findstat(41).set_description("")                              # optional -- internet
            sage: findstat(41).modified()                                       # optional -- internet
            True

        """
        return self._modified

    def collection(self):
        r"""
        Return the FindStat collection of the statistic.

        OUTPUT:

        The FindStat collection of the statistic as an instance of
        :class:`FindStatCollection`.

        EXAMPLES::

            sage: findstat(1).collection()                                      # optional -- internet
            Cc0001: Permutations
        """
        return self._collection

    def function(self):
        r"""
        Return the function used to compute the values of the statistic.

        OUTPUT:

        The function used to compute the values of the statistic, or
        ``None``.

        EXAMPLES::

            sage: findstat("Permutations", lambda pi: pi.length()).function()   # optional -- internet
            ...
            <function <lambda> at ...>
        """
        return self._function

    def first_terms(self):
        r"""
        Return the first terms of the statistic.

        OUTPUT:

        A list of pairs of the form ``(object, value)`` where
        ``object`` is a sage object representing an element of the
        appropriate collection and ``value`` is an integer.  If the
        statistic is in the FindStat database, the list contains
        exactly the pairs in the database.

        EXAMPLES::

            sage: findstat(1).first_terms()                                     # optional -- internet,random
            [([1], 1),
             ([1, 2], 1),
             ([2, 1], 1),
             ([1, 2, 3], 1),
             ([1, 3, 2], 1),
             ([2, 1, 3], 1),
             ...

        TESTS::

            sage: r = findstat({d: randint(1,1000) for d in DyckWords(4)}); r   # optional -- internet
            a new statistic on Cc0005: Dyck paths

            sage: isinstance(r.first_terms(), list)                             # optional -- internet
            True
            sage: all(isinstance(e, tuple) and len(e)==2 and isinstance(e[1], (ZZ, Integer, int)) for e in r.first_terms())      # optional -- internet
            True
        """
        return self._first_terms

    def first_terms_str(self):
        r"""
        Return the first terms of the statistic in the format needed
        for a FindStat query.

        OUTPUT:

        A string, where each line is of the form ``object => value``,
        where ``object`` is the string representation of an element
        of the appropriate collection as used by FindStat and value
        is an integer.

        EXAMPLES::

            sage: findstat(1).first_terms_str()[:10]                            # optional -- internet,random
            '[1] => 1\r\n'

        """
        if self._first_terms is not None:
            to_str = self._collection.to_string()
            return "\r\n".join([to_str(key) + " => " + str(val)
                                for (key, val) in self._first_terms])
        return ""

    # apart from efficiency considerations, it is important to make
    # this lazy because the method _get_level is not available for
    # unsupported collections.
    @lazy_attribute
    def _generating_functions_dict(self):
        r""" Return the generating functions of ``self`` in a dictionary,
        computed from ``self._data``.

        TESTS:

            sage: s = findstat((DyckWords(4), range(14))); s                    # optional -- internet, indirect doctest
            a new statistic on Cc0005: Dyck paths
            sage: s.generating_functions()                                      # optional -- internet, indirect doctest
            {4: q^13 + q^12 + q^11 + q^10 + q^9 + q^8 + q^7 + q^6 + q^5 + q^4 + q^3 + q^2 + q + 1}

            sage: C = AlternatingSignMatrices(4)
            sage: s = findstat((C, range(C.cardinality()))); s                  # optional -- internet, indirect doctest
            a new statistic on Cc0017: Alternating sign matrices
            sage: s.generating_functions()                                      # optional -- internet, indirect doctest
            {4: q^41 + q^40 + q^39 + q^38 + q^37 + q^36 + q^35 + q^34 + q^33 + q^32 + q^31 + q^30 + q^29 + q^28 + q^27 + q^26 + q^25 + q^24 + q^23 + q^22 + q^21 + q^20 + q^19 + q^18 + q^17 + q^16 + q^15 + q^14 + q^13 + q^12 + q^11 + q^10 + q^9 + q^8 + q^7 + q^6 + q^5 + q^4 + q^3 + q^2 + q + 1}
        """
        c = self.collection()
        gfs = dict()
        for (elements, elements_str, values) in self._data:
            l = c._get_level(elements[0])
            if l in c._levels and all(l == c._get_level(e) for e in elements[1:]):
                gfs[l] = gfs.get(l, dict())
                for v in values:
                    gfs[l][v] = gfs[l].get(v, 0) + 1

        for (l, p) in gfs.items():
            ns = sum(p.values())
            nc = c._levels[l]
            if ns < nc:
                del gfs[l]
            elif ns > nc:
                raise ValueError("FindStat delivers more statistic values than the level has elements.  This should not happen.  Please send an email to the developers.")

        return gfs

    def generating_functions(self, style="polynomial"):
        r"""
        Return the generating functions of ``self`` in a dictionary.

        The keys of this dictionary are the levels for which the
        generating function of ``self`` can be computed from the data
        of this statistic, and each value represents a generating
        function for one level, as a polynomial, as a dictionary, or as
        a list of coefficients.

        INPUT:

        - a string -- (default:"polynomial") can be
          "polynomial", "dictionary", or "list".

        OUTPUT:

        - if ``style`` is ``"polynomial"``, the generating function is
          returned as a polynomial.

        - if ``style`` is ``"dictionary"``, the generating function is
          returned as a dictionary representing the monomials of the
          generating function.

        - if ``style`` is ``"list"``, the generating function is
          returned as a list of coefficients of the generating function.

        EXAMPLES::

            sage: st = findstat(18)                                             # optional -- internet

            sage: st.generating_functions()                                     # optional -- internet,random
            {2: q + 1,
             3: q^3 + 2*q^2 + 2*q + 1,
             4: q^6 + 3*q^5 + 5*q^4 + 6*q^3 + 5*q^2 + 3*q + 1,
             5: q^10 + 4*q^9 + 9*q^8 + 15*q^7 + 20*q^6 + 22*q^5 + 20*q^4 + 15*q^3 + 9*q^2 + 4*q + 1,
             6: q^15 + 5*q^14 + 14*q^13 + 29*q^12 + 49*q^11 + 71*q^10 + 90*q^9 + 101*q^8 + 101*q^7 + 90*q^6 + 71*q^5 + 49*q^4 + 29*q^3 + 14*q^2 + 5*q + 1}

            sage: st.generating_functions(style="dictionary")                   # optional -- internet,random
            {2: {0: 1, 1: 1},
             3: {0: 1, 1: 2, 2: 2, 3: 1},
             4: {0: 1, 1: 3, 2: 5, 3: 6, 4: 5, 5: 3, 6: 1},
             5: {0: 1, 1: 4, 2: 9, 3: 15, 4: 20, 5: 22, 6: 20, 7: 15, 8: 9, 9: 4, 10: 1},
             6: {0: 1,
              1: 5,
              2: 14,
              3: 29,
              4: 49,
              5: 71,
              6: 90,
              7: 101,
              8: 101,
              9: 90,
              10: 71,
              11: 49,
              12: 29,
              13: 14,
              14: 5,
              15: 1}}

            sage: st.generating_functions(style="list")                         # optional -- internet,random
            {2: [1, 1],
             3: [1, 2, 2, 1],
             4: [1, 3, 5, 6, 5, 3, 1],
             5: [1, 4, 9, 15, 20, 22, 20, 15, 9, 4, 1],
             6: [1, 5, 14, 29, 49, 71, 90, 101, 101, 90, 71, 49, 29, 14, 5, 1]}
        """
        gfs = self._generating_functions_dict
        if style == "dictionary":
            return gfs
        elif style == "list":
            return { key : [ gfs[key][deg] if deg in gfs[key] else 0
                             for deg in range(min(gfs[key]),max(gfs[key])+1) ]
                     for key in gfs }
        elif style == "polynomial":
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            from sage.rings.integer_ring import ZZ
            P = PolynomialRing(ZZ,"q")
            q = P.gen()
            return { level : sum( coefficient * q**exponent
                                  for exponent,coefficient in gen_dict.iteritems() )
                     for level, gen_dict in gfs.iteritems()}
        else:
            raise ValueError("The argument 'style' (='%s') must be 'dictionary', 'polynomial', or 'list'."%style)

    def oeis_search(self, search_size=32, verbose=True):
        r"""
        Search the OEIS for the generating function of the statistic.

        INPUT:

        - ``search_size`` (default:32) the number of integers in the
          sequence. If too big, the OEIS result is corrupted.

        - ``verbose`` (default:True) if true, some information about
          the search are printed.

        OUTPUT:

        - a tuple of OEIS sequences, see
          :meth:`sage.databases.oeis.OEIS.find_by_description` for more
          information.

        EXAMPLES::

            sage: st = findstat(18)                                             # optional -- internet

            sage: st.oeis_search()                                              # optional -- internet,random
            Searching the OEIS for "1,1  1,2,2,1  1,3,5,6,5,3,1  1,4,9,15,20,22,20,15,9,4,1  1,5,14,29,49,71,90,101"

            0: A008302: Triangle of Mahonian numbers T(n,k): coefficients in expansion of Product_{i=0..n-1} (1 + x + ... + x^i), where k ranges from 0 to A000217(n-1).

            sage: st.oeis_search(search_size=13)                                # optional -- internet,random
            Searching the OEIS for "1,1  1,2,2,1  1,3,5,6,5,3,1"

            0: A008302: Triangle of Mahonian numbers T(n,k): coefficients in expansion of Product_{i=0..n-1} (1 + x + ... + x^i), where k ranges from 0 to A000217(n-1).
            1: A115570: Array read by rows: row n (n>= 1) gives the Betti numbers for the n-th element of the Weyl group of type A3 (in Goresky's standard ordering).
            2: A187447: Array for all multiset choices (multiset repetition class representatives in Abramowitz-Stegun order).
        """
        from sage.databases.oeis import oeis

        gen_funcs = self.generating_functions(style="list")

        OEIS_string = ""
        keys = sorted(gen_funcs.keys())
        counter = 0
        for key in keys:
            gen_func = gen_funcs[key]
            while gen_func[0] == 0:
                gen_func.pop(0)
            # we strip the result according to the search size. -- stumpc5, 2015-09-27
            gen_func = gen_func[:search_size]
            counter += len(gen_func)
            if search_size > 0:
                search_size -= len(gen_func)
            OEIS_func_string     = ",".join( str(coefficient) for coefficient in gen_func )
            OEIS_string         += OEIS_func_string + "  "
        OEIS_string = OEIS_string.strip()
        if counter >= 4:
            if verbose:
                print('Searching the OEIS for "%s"' % OEIS_string)
            return oeis( OEIS_string )
        else:
            if verbose:
                print("Too little information to search the OEIS for this statistic (only %s values given)." % counter)
            return

    def description(self):
        r"""
        Return the description of the statistic.

        OUTPUT:

        A string, whose first line is used as the name of the
        statistic.

        EXAMPLES::

            sage: print(findstat(1).description())    # optional -- internet,random
            The number of ways to write a permutation as a minimal length product of simple transpositions.
            <BLANKLINE>
            That is, the number of reduced words for the permutation.  E.g., there are two reduced words for $[3,2,1] = (1,2)(2,3)(1,2) = (2,3)(1,2)(2,3)$.
        """
        return self._description

    def set_description(self, value):
        r"""
        Set the description of the statistic.

        INPUT:

        - a string -- the name of the statistic followed by its
          description on a separate line.

        OUTPUT:

        - Raise an error, if the query has a match with no
          intermediate combinatorial maps.

        This information is used when submitting the statistic with
        :meth:`submit`.

        EXAMPLES::

            sage: s = findstat([(d, randint(1,1000)) for d in DyckWords(4)]); s # optional -- internet
            a new statistic on Cc0005: Dyck paths
            sage: s.set_description("Random values on Dyck paths.\r\nNot for submission.") # optional -- internet
            sage: s                                                             # optional -- internet
            a new statistic on Cc0005: Dyck paths
            sage: s.name()                                                      # optional -- internet
            'Random values on Dyck paths.'
            sage: print(s.description())                              # optional -- internet
            Random values on Dyck paths.
            Not for submission.
        """
        self._raise_error_modifying_statistic_with_perfect_match()

        if value != self._description:
            self._modified = True
            self._description = value

    def name(self):
        r"""
        Return the name of the statistic.

        OUTPUT:

        A string, which is just the first line of the description of
        the statistic.

        EXAMPLES::

            sage: findstat(1).name()                                            # optional -- internet,random
            u'The number of ways to write a permutation as a minimal length product of simple transpositions.'
        """
        # this needs to be decided how to do properly
        if hasattr(self,"_name"):
            return self._name
        else:
            return self._description.partition(FINDSTAT_SEPARATOR_NAME)[0]

    def references(self):
        r"""
        Return the references associated with the statistic.

        OUTPUT:

        An instance of :class:`sage.databases.oeis.FancyTuple`, each
        item corresponds to a reference.

        .. TODO::

            Since the references in the database are sometimes not
            formatted properly, this method is unreliable.  The
            string representation can be obtained via
            :attr:`_references`.

        EXAMPLES::

            sage: findstat(1).references()                                      # optional -- internet,random
            0: P. Edelman and C. Greene, Balanced tableaux, Adv. in Math., 63 (1987), pp. 42-99.
            1: [[OEIS:A005118]]
            2: [[oeis:A246865]]
        """
        l = [ref.strip() for ref in self._references.split(FINDSTAT_SEPARATOR_REFERENCES)]
        return FancyTuple([ref for ref in l if ref != ""])

    def set_references(self, value):
        r"""
        Set the references associated with the statistic.

        INPUT:

        - a string -- the individual references should be separated
          by FINDSTAT_SEPARATOR_REFERENCES, which is "\\r\\n".

        OUTPUT:

        - Raise an error, if the query has a match with no
          intermediate combinatorial maps.

        This information is used when submitting the statistic with
        :meth:`submit`.

        EXAMPLES::

            sage: s = findstat([(d, randint(1,1000)) for d in DyckWords(4)]); s # optional -- internet
            a new statistic on Cc0005: Dyck paths
            sage: s.set_references("[1] The wonders of random Dyck paths, Anonymous Coward, [[arXiv:1102.4226]].\r\n[2] [[oeis:A000001]]")  # optional -- internet
            sage: s.references()                                                # optional -- internet
            0: [1] The wonders of random Dyck paths, Anonymous Coward, [[arXiv:1102.4226]].
            1: [2] [[oeis:A000001]]

        """
        self._raise_error_modifying_statistic_with_perfect_match()

        if value != self._references:
            self._modified = True
            self._references = value

    def code(self):
        r"""
        Return the code associated with the statistic.

        OUTPUT:

        A string.  Contributors are encouraged to submit sage code in the form::

            def statistic(x):
                ...

        but the string may also contain code for other computer
        algebra systems.

        EXAMPLES::

            sage: print(findstat(1).code())                                   # optional -- internet,random
            def statistic(x):
                return len(x.reduced_words())

            sage: print(findstat(118).code())                                 # optional -- internet,random
            (* in Mathematica *)
            tree = {{{{}, {}}, {{}, {}}}, {{{}, {}}, {{}, {}}}};
            Count[tree, {{___}, {{___}, {{___}, {___}}}}, {0, Infinity}]
        """
        return self._code

    def set_code(self, value):
        r"""
        Set the code associated with the statistic.

        INPUT:

        - a string -- contributors are encouraged to submit sage code
          in the form::

            def statistic(x):
                ...

        OUTPUT:

        - Raise an error if the query has a match with no
          intermediate combinatorial maps.

        This information is used when submitting the statistic with
        :meth:`submit`.

        EXAMPLES::

            sage: s = findstat([(d, randint(1,1000)) for d in DyckWords(4)])    # optional -- internet
            sage: s.set_code("def statistic(x):\r\n    return randint(1,1000)") # optional -- internet
            sage: print(s.code())                 # optional -- internet
            def statistic(x):
                return randint(1,1000)
        """
        self._raise_error_modifying_statistic_with_perfect_match()

        if value != self._code:
            self._modified = True
            self._code = value

    ######################################################################
    # browse current statistic
    ######################################################################

    def browse(self):
        r"""
        Open the FindStat web page of the statistic in a browser.

        EXAMPLES::

            sage: findstat(41).browse()                                         # optional -- webbrowser
        """
        if self._query == "ID":
            webbrowser.open(FINDSTAT_URL_BROWSE + self.id_str())
        else:
            raise NotImplementedError("It would be nice to show the result of a FindStat query in the webbrowser, but we do not know how to do this yet.")

    ######################################################################
    # submit current (possibly incompletely defined) statistic
    ######################################################################

    def submit(self, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Open the FindStat web page for editing the statistic in a browser.

        INPUT:

        - ``max_values`` -- integer (default:
          ``FINDSTAT_MAX_SUBMISSION_VALUES``); if :meth:`function` is
          defined and the statistic is a new statistic, use
          :meth:`FindStatCollection.first_terms` to produce at most
          ``max_values`` terms.

        OUTPUT:

        - Raise an error if the query has a match with no
          intermediate combinatorial maps.

        EXAMPLES::

            sage: s = findstat(DyckWords(4), lambda x: randint(1,1000)); s      # optional -- internet
            a new statistic on Cc0005: Dyck paths

        The following uses ``lambda x: randint(1,1000)`` to produce
        14 terms, because ``min(DyckWords(4).cardinality(),
        FINDSTAT_MAX_SUBMISSION_VALUES)`` is 14::

            sage: s.submit()                                                    # optional -- webbrowser

        """
        self._raise_error_modifying_statistic_with_perfect_match()

        # if the statistic is given as a function, and we have a new
        # statistic then update first_terms

        # it is not clear whether we want to do this also for old statistics.
        if self.function() and self.id() == 0:
            self._first_terms = self.collection().first_terms(self.function(),
                                                              max_values=max_values)

        args = dict()
        args[FINDSTAT_STATISTIC_IDENTIFIER]  = self._id
        args[FINDSTAT_STATISTIC_COLLECTION]  = str(self.collection().id())
        args[FINDSTAT_STATISTIC_DATA]        = self.first_terms_str()
        args[FINDSTAT_STATISTIC_DESCRIPTION] = self._description
        args[FINDSTAT_STATISTIC_REFERENCES]  = self._references
        args[FINDSTAT_STATISTIC_CODE]        = self.code()
        args[FINDSTAT_POST_SAGE_CELL]        = ""
        args[FINDSTAT_POST_EDIT]             = ""
        args[FINDSTAT_POST_AUTHOR]           = findstat._user_name
        args[FINDSTAT_POST_EMAIL]            = findstat._user_email

        assert set(args.keys()) == FINDSTAT_EDIT_FIELDS, "It appears that the list of required post variables for editing a statistic has changed.  Please update FindStatStatistic.submit()."

        # write the file
        f = tempfile.NamedTemporaryFile(delete=False)
        _ = verbose("Created temporary file %s" %f.name, caller_name='FindStat')
        f.write(FINDSTAT_POST_HEADER)
        if self.id() == 0:
            f.write(FINDSTAT_NEWSTATISTIC_FORM_HEADER %FINDSTAT_URL_NEW)
        else:
            f.write(FINDSTAT_NEWSTATISTIC_FORM_HEADER %(FINDSTAT_URL_EDIT+self.id_str()))
        for key, value in args.iteritems():
            _ = verbose("writing argument %s" %key, caller_name='FindStat')
            value_encoded = cgi.escape(str(value), quote=True)
            _ = verbose("%s" %value_encoded, caller_name='FindStat')
            f.write((FINDSTAT_NEWSTATISTIC_FORM_FORMAT %(key, value_encoded)))
        f.write(FINDSTAT_NEWSTATISTIC_FORM_FOOTER)
        f.close()
        _ = verbose("Opening file with webbrowser", caller_name='FindStat')
        _ = webbrowser.open(f.name)

        _ = verbose("Waiting a little before deleting the temporary file", caller_name='FindStat')
        time.sleep(1)

        f.unlink(f.name)

    # editing and submitting is really the same thing
    edit = submit

# helper for generation of CartanTypes
def _finite_irreducible_cartan_types_by_rank(n):
    """
    Return the Cartan types of rank n.

    INPUT:

    - n -- an integer.

    OUTPUT:

    The list of Cartan types of rank n.

    TESTS::

        sage: from sage.databases.findstat import _finite_irreducible_cartan_types_by_rank
        sage: _finite_irreducible_cartan_types_by_rank(2)
        [['A', 2], ['B', 2], ['G', 2]]
    """
    cartan_types      = [ CartanType(['A',n]) ]
    if n >= 2:
        cartan_types += [ CartanType(['B',n]) ]
    if n >= 3:
        cartan_types += [ CartanType(['C',n]) ]
    if n >= 4:
        cartan_types += [ CartanType(['D',n]) ]
    if n in [6,7,8]:
        cartan_types += [ CartanType(['E',n]) ]
    if n == 4:
        cartan_types += [ CartanType(['F',n]) ]
    if n == 2:
        cartan_types += [ CartanType(['G',n]) ]
    return cartan_types

class FindStatCollection(Element):
    r"""
    A FindStat collection.

    :class:`FindStatCollection` is a class representing a
    combinatorial collection available in the FindStat database.

    Its main use is to allow easy specification of the combinatorial
    collection when using :class:`findstat<FindStat>`.  It also
    provides methods to quickly access its FindStat web page
    (:meth:`browse`), check whether a particular element is actually
    in the range considered by FindStat (:meth:`in_range`), etc.

    INPUT:

    One of the following:

    - a string eg. 'Dyck paths' or 'DyckPaths', case-insensitive, or

    - an integer designating the FindStat id of the collection, or

    - a sage object belonging to a collection, or

    - an iterable producing a sage object belonging to a collection.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatCollection
        sage: FindStatCollection("Dyck paths")                                  # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(5)                                             # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(DyckWord([1,0,1,0]))                           # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(DyckWords(2))                                  # optional -- internet
        Cc0005: Dyck paths

    SEEALSO:

    :class:`FindStatCollections`

    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, entry):
        """
        Retrieve a collection from the database.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection(0)                                         # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: Could not find FindStat collection for 0.
        """
        return FindStatCollections()(entry)

    def __init__(self, parent, id, c, sageconstructor_overridden):
        """Initialize the collection.

        This should only be called in
        :meth:`FindStatCollections()._element_constructor_` via
        `element_class`.

        INPUT:

        - ``parent`` -- :class:`FindStatCollections`.

        - ``id`` -- the FindStat identifier of the collection.

        - ``c`` -- a tuple containing the properties of the
          collection, such as its name, the corresponding class in
          sage, and so on.

        - ``sageconstructor_overridden`` -- either ``None`` or an
          iterable which yields a subset of the elements of the
          collection.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection(5).parent()                                # optional -- internet
            Set of combinatorial collections used by FindStat

        """
        self._id = id
        (self._name, self._name_plural, self._url_name,
         self._sageclass, self._sageconstructor, self._levels, self._get_level,
         self._in_range, self._to_str, self._from_str) = c
        self._sageconstructor_overridden = sageconstructor_overridden

        Element.__init__(self, parent)

    def __reduce__(self):
        """Return a function and its arguments needed to create this
        collection.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("Permutations")                        # optional -- internet
            sage: loads(dumps(c)) == c                                          # optional -- internet
            True

        """
        return (FindStatCollection, (self.id(),))

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatCollection, FindStatCollections
            sage: FindStatCollection("Permutations") == FindStatCollection("Permutations")          # optional -- internet
            True

            sage: FindStatCollection("Permutations") == FindStatCollection("Integer Partitions")    # optional -- internet
            False

            sage: FindStatCollection("Permutations") != FindStatCollection("Permutations")          # optional -- internet
            False

            sage: FindStatCollection("Permutations") != FindStatCollection("Integer Partitions")    # optional -- internet
            True

            sage: FindStatCollection("Permutations") == 1                                           # optional -- internet
            False

            sage: FindStatCollection("Permutations") != 1                                           # optional -- internet
            True

            sage: sorted(c for c in FindStatCollections())[0]                                       # optional -- internet
            Cc0001: Permutations
        """
        return self.id().__cmp__(other.id())

    def is_supported(self):
        """
        Check whether the collection is fully supported by the interface.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection(1).is_supported()                          # optional -- internet
            True

            sage: FindStatCollection(24).is_supported()                         # optional -- internet, random
            False

        """
        try:
            self._sageconstructor(self._levels.keys()[0])
            return True
        except NotImplementedError:
            return False

    def in_range(self, element):
        r"""
        Check whether an element of the collection is in FindStat's precomputed range.

        INPUT:

        - ``element`` -- a sage object that belongs to the collection.

        OUTPUT:

        ``True``, if ``element`` is used by the FindStat search
        engine, and ``False`` if it is ignored.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.in_range(GelfandTsetlinPattern([[2, 1], [1]]))              # optional -- internet
            True
            sage: c.in_range(GelfandTsetlinPattern([[3, 1], [1]]))              # optional -- internet
            True
            sage: c.in_range(GelfandTsetlinPattern([[4, 1], [1]]))              # optional -- internet,random
            False


        TESTS::

            sage: from sage.databases.findstat import FindStatCollections
            sage: l = FindStatCollections()                                     # optional -- internet
            sage: long = [9, 12, 14, 20]
            sage: for c in l:                                                   # optional -- internet, random
            ....:     if c.id() not in long and c.is_supported():
            ....:         f = c.first_terms(lambda x: 1, max_values=10000)
            ....:         print("{} {} {}".format(c, len(f), all(c.in_range(e) for e, _ in f)))
            ....:
            Cc0001: Permutations 10000 True
            Cc0002: Integer partitions 270 True
            Cc0005: Dyck paths 2054 True
            Cc0006: Integer compositions 510 True
            Cc0007: Standard tableaux 3734 True
            Cc0010: Binary trees 2054 True
            Cc0013: Cores 100 True
            Cc0017: Alternating sign matrices 7916 True
            Cc0018: Gelfand-Tsetlin patterns 934 True
            Cc0019: Semistandard tableaux 10000 True
            Cc0021: Ordered trees 2055 True
            Cc0022: Finite Cartan types 31 True
            Cc0023: Parking functions 10000 True

        """
        return self._in_range(element, self._levels.keys())

    def first_terms(self, statistic, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Compute the first few terms of the given statistic.

        INPUT:

        - ``statistic`` -- a callable.

        - ``max_values`` -- the number of terms to compute at most.

        OUTPUT:

        A list of pairs of the form (object, value).

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.first_terms(lambda x: 1, max_values=10)                     # optional -- internet,random
            [([[0]], 1),
             ([[1]], 1),
             ([[2]], 1),
             ([[3]], 1),
             ([[0, 0], [0]], 1),
             ([[1, 0], [0]], 1),
             ([[1, 0], [1]], 1),
             ([[1, 1], [1]], 1),
             ([[2, 0], [0]], 1),
             ([[2, 0], [1]], 1)]
        """
        if self._sageconstructor_overridden is None:
            g = (x for n in self._levels.keys() for x in self._sageconstructor(n))
        else:
            g = self._sageconstructor_overridden

        return [(x, statistic(x)) for (x,_) in zip(g, xrange(max_values))]

    def id(self):
        r"""
        Return the FindStat identifier of the collection.

        OUTPUT:

        The FindStat identifier of the collection as an integer.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.id()                                                        # optional -- internet
            18
        """
        return self._id

    def id_str(self):
        r"""
        Return the FindStat identifier of the collection.

        OUTPUT:

        The FindStat identifier of the collection as a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.id_str()                                                    # optional -- internet
            'Cc0018'
        """
        id = str(self.id())
        return 'Cc0000'[:-len(id)] + id

    def browse(self):
        r"""
        Open the FindStat web page of the collection in a browser.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Permutations").browse()                   # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL + self._url_name)

    def to_string(self):
        r"""
        Return a function that returns the FindStat normal
        representation given an object.

        OUTPUT:

        The function that produces the string representation as
        needed by the FindStat search webpage.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: p = Poset((range(3), [[0, 1], [1, 2]]))                       # optional -- internet
            sage: c = FindStatCollection("Posets")                              # optional -- internet
            sage: c.to_string()(p)                                              # optional -- internet
            '([(0, 2), (2, 1)], 3)'
        """
        return self._to_str

    def from_string(self):
        r"""
        Return a function that returns the object given the FindStat
        normal representation.

        OUTPUT:

        The function that produces the sage object given its FindStat
        normal representation as a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("Posets")                              # optional -- internet
            sage: p = c.from_string()('([(0, 2), (2, 1)], 3)')                  # optional -- internet
            sage: p.cover_relations()                                           # optional -- internet
            [[0, 2], [2, 1]]

            sage: c = FindStatCollection("Binary Words")                        # optional -- internet
            sage: w = c.from_string()('010101')                                 # optional -- internet
            sage: w in c._sageconstructor(6)                                    # optional -- internet
            True
        """
        return self._from_str

    def _repr_(self):
        r"""
        Return the representation of the FindStat collection.

        OUTPUT:

        The representation, including the identifier and the name.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Binary trees")                            # optional -- internet
            Cc0010: Binary trees
        """
        return "%s: %s" %(self.id_str(), self._name_plural)

    def name(self, style="singular"):
        r"""
        Return the name of the FindStat collection.

        INPUT:

        - a string -- (default:"singular") can be
          "singular", or "plural".

        OUTPUT:

        The name of the FindStat collection, in singular or in plural.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Binary trees").name()                     # optional -- internet
            u'Binary tree'

            sage: FindStatCollection("Binary trees").name(style="plural")       # optional -- internet
            u'Binary trees'
        """
        if style == "singular":
            return self._name
        elif style == "plural":
            return self._name_plural
        else:
            raise ValueError("Argument 'style' (=%s) must be 'singular' or 'plural'."%style)

class FindStatCollections(Parent, UniqueRepresentation):
    r"""
    The class of FindStat collections.

    The elements of this class are combinatorial collections in
    FindStat as of August 2015.  If a new collection was added to the
    web service since then, the dictionary ``_findstat_collections``
    in this class has to be updated accordingly.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatCollections
        sage: sorted(c for c in FindStatCollections())                          # optional -- internet, random
        [Cc0001: Permutations,
         Cc0002: Integer partitions,
         Cc0005: Dyck paths,
         Cc0006: Integer compositions,
         Cc0007: Standard tableaux,
         Cc0009: Set partitions,
         Cc0010: Binary trees,
         Cc0012: Perfect matchings,
         Cc0013: Cores,
         Cc0014: Posets,
         Cc0017: Alternating sign matrices,
         Cc0018: Gelfand-Tsetlin patterns,
         Cc0019: Semistandard tableaux,
         Cc0020: Graphs,
         Cc0021: Ordered trees,
         Cc0022: Finite Cartan types,
         Cc0023: Parking functions]
    """

    # we set up a dict of FindStat collections containing, with key
    # being the FINDSTAT_COLLECTION_IDENTIFIER to tuples, containing in this order:

    # * the FindStat name                  (None)                (FINDSTAT_COLLECTION_NAME)
    # * the FindStat name plural           (None)                (FINDSTAT_COLLECTION_NAME_PLURAL)
    # * url's as needed by FindStat        (None)                (FINDSTAT_COLLECTION_NAME_WIKI)
    # * sage element constructor (used in _element_constructor_ to check whether an element is an instance of this collection, this should be replaced by a function that recognises elements)
    # * sage constructor                                         (would be parent_initializer, accepts a level to return the component of this degree)
    # * a dictionary from levels to sizes  (None)                (FINDSTAT_COLLECTION_LEVELS)
    #   the levels are used as arguments for the sage constructor
    # * a function to check whether an object is produced by applying the constructor to some element in the list
    # * the (FindStat) string representations of the sage object (would be element_repr)
    # * sage constructors of the FindStat string representation  (would be element_constructor)

    # when adding a collection, note the following:

    # to recognize a collection given a sage object,
    # :meth:`_element_constructor_` checks whether the object is an
    # instance of one of the sage element constructors.

    # sage objects are normalized using the method :meth:`to_string`.
    # This method should apply to objects produced by
    # :meth:`first_terms` as well as to objects produced by
    # :meth:`from_string`.

    # several fields are initialised with 'None', they are updated
    # upon the first call to this class
    _findstat_collections = {
        24: [None, None, None, Word_class, lambda x: Words([0,1], x), None,
             lambda x: x.length(),
             lambda x, l: x.length() in l,
             str,
             lambda x: Word((int(e) for e in str(x)), alphabet=[0,1])],
        17: [None, None, None, AlternatingSignMatrix, AlternatingSignMatrices, None,
             lambda x: x.to_matrix().nrows(),
             lambda x, l: x.to_matrix().nrows() in l,
             lambda x: str(map(list, list(x._matrix))),
             lambda x: AlternatingSignMatrix(literal_eval(x))],
        10: [None, None, None, BinaryTree,            BinaryTrees,             None,
             lambda x: x.node_number(),
             lambda x, l: x.node_number() in l,
             str,
             lambda x: BinaryTree(str(x))],
        13: [None, None, None, Core,                  lambda x: Cores(x[1], x[0]),
             None,
             lambda x: (x.length(), x.k()),
             lambda x, l: (x.length(), x.k()) in l,
             lambda X: "( "+X._repr_()+", "+str(X.k())+" )",
             lambda x: (lambda pi, k: Core(pi, k))(*literal_eval(x))],
        5:  [None, None, None, DyckWord,              DyckWords,               None,
             lambda x: (x.length()/2),
             lambda x, l: (x.length()/2) in l,
             lambda x: str(list(DyckWord(x))),
             lambda x: DyckWord(literal_eval(x))],
        22: [None, None, None, CartanType_abstract,   _finite_irreducible_cartan_types_by_rank,
             None,
             lambda x: x.rank(),
             lambda x, l: x.rank() in l,
             str,
             lambda x: CartanType(*literal_eval(str(x)))],
        18: [None, None, None, GelfandTsetlinPattern, lambda x: GelfandTsetlinPatterns(*x),
             None,
             None,
             lambda x, l: any(len(x) == s and max([0] + [max(row) for row in x]) <= m for s, m in l),
             str,
             lambda x: GelfandTsetlinPattern(literal_eval(x))],
        20: [None, None, None, Graph,                 graphs,
             None,
             lambda x: x.num_verts(),
             lambda x, l: x.num_verts() in l,
             lambda X: str((sorted(X.canonical_label().edges(False)), X.num_verts())),
             lambda x: (lambda E, V: Graph([range(V), lambda i,j: (i,j) in E or (j,i) in E], immutable=True))(*literal_eval(x))],
        6:  [None, None, None, Composition,           Compositions,            None,
             lambda x: x.size(),
             lambda x, l: x.size() in l,
             str,
             lambda x: Composition(literal_eval(x))],
        2:  [None, None, None, Partition,             Partitions,              None,
             lambda x: x.size(),
             lambda x, l: x.size() in l,
             str,
             lambda x: Partition(literal_eval(x))],
        21: [None, None, None, OrderedTree,           OrderedTrees,            None,
             lambda x: x.node_number(),
             lambda x, l: x.node_number() in l,
             str,
             lambda x: OrderedTree(literal_eval(x))],
        23: [None, None, None, ParkingFunction_class, ParkingFunctions,        None,
             lambda x: len(x),
             lambda x, l: len(x) in l,
             str,
             lambda x: ParkingFunction(literal_eval(x))],
        12: [None, None, None, PerfectMatching,       PerfectMatchings,        None,
             lambda x: x.size(),
             lambda x, l: x.size() in l,
             str,
             lambda x: PerfectMatching(literal_eval(x))],
        1:  [None, None, None, Permutation,           Permutations,            None,
             lambda x: x.size(),
             lambda x, l: x.size() in l,
             str,
             lambda x: Permutation(literal_eval(x))],
        14: [None, None, None, FinitePoset,           posets,                  None,
             lambda x: x.cardinality(),
             lambda x, l: x.cardinality() in l,
             lambda X: str((sorted(X._hasse_diagram.canonical_label().cover_relations()), len(X._hasse_diagram.vertices()))),
             lambda x: (lambda R, E: Poset((range(E), R)))(*literal_eval(x))],
        19: [None, None, None, SemistandardTableau,   lambda x: SemistandardTableaux(x),
             None,
             lambda x: x.size(),
             lambda x, l: x.size() in l,
             str,
             lambda x: SemistandardTableau(literal_eval(x))],
        9:  [None, None, None, SetPartition,          SetPartitions,           None,
             lambda x: x.size(),
             lambda x, l: x.size() in l,
             str,
             lambda x: SetPartition(literal_eval(x.replace('{','[').replace('}',']')))],
        7:  [None, None, None, StandardTableau,       StandardTableaux,        None,
             lambda x: x.size(),
             lambda x, l: x.size() in l,
             str,
             lambda x: StandardTableau(literal_eval(x))]}

    def _raise_unsupported_error(*args):
        """A placeholder function for unsupported collections that raises an
        error.

        sage: from sage.databases.findstat import FindStatCollection
        sage: FindStatCollection(24).first_terms(lambda x: 1);                  # optional -- internet,random,indirect doctest
        Traceback (most recent call last):
        ...
        NotImplementedError: This FindStatCollection is not yet supported.
        """
        raise NotImplementedError("This FindStatCollection is not yet supported.")

    # The following is used for unknown collections, with the
    # intention to make as much as possible working.
    #
    # In particular, we keep the elements provided by findstat.org as
    # strings.  Therefore, to_string and from_string should do
    # nothing.  Furthermore, in_range always returns True to make
    # submission of sequences possible.
    _findstat_unsupported_collection_default = [None, None, None,
                                                _raise_unsupported_error,
                                                _raise_unsupported_error,
                                                None,
                                                _raise_unsupported_error,
                                                lambda x, l: True,
                                                lambda x: x, lambda x: x]

    def __init__(self):
        """
        Fetch the collections from FindStat.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollections
            sage: C = FindStatCollections()                                     # optional -- internet
            sage: TestSuite(C).run()                                            # optional -- internet
        """
        for j in json.load(urlopen(FINDSTAT_URL_DOWNLOADS_COLLECTIONS)):
            id = j[FINDSTAT_COLLECTION_IDENTIFIER]
            if id in self._findstat_collections:
                c = self._findstat_collections[id]
            else:
                print("There is a new collection available at `%s`: %s." % (findstat, j[FINDSTAT_COLLECTION_NAME_PLURAL]))
                print("To use it with this interface, it has to be added to the dictionary FindStatCollections._findstat_collections in src/sage/databases/findstat.py of the SageMath distribution.  Please open a ticket on trac!")
                self._findstat_collections[id] = list(self._findstat_unsupported_collection_default)
                c = self._findstat_collections[id]
            c[0] = j[FINDSTAT_COLLECTION_NAME]
            c[1] = j[FINDSTAT_COLLECTION_NAME_PLURAL]
            c[2] = j[FINDSTAT_COLLECTION_NAME_WIKI]
            c[5] =  {literal_eval(key):value for key,value in
                     j[FINDSTAT_COLLECTION_LEVELS].iteritems()}

        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, entry):
        """Retrieve a FindStat collection from the database.

        INPUT:

        see :class:`FindStatCollection`.

        TESTS:

        Create an object and find its collection::

            sage: from sage.databases.findstat import FindStatCollection, FindStatCollections
            sage: [FindStatCollection(c.first_terms(lambda x: 0, max_values=1)[0][0]) for c in FindStatCollections() if c.is_supported()]   # optional -- internet, random
            [Cc0001: Permutations,
             Cc0002: Integer partitions,
             Cc0005: Dyck paths,
             Cc0006: Integer compositions,
             Cc0007: Standard tableaux,
             Cc0009: Set partitions,
             Cc0010: Binary trees,
             Cc0012: Perfect matchings,
             Cc0013: Cores,
             Cc0014: Posets,
             Cc0017: Alternating sign matrices,
             Cc0018: Gelfand-Tsetlin patterns,
             Cc0019: Semistandard tableaux,
             Cc0020: Graphs,
             Cc0021: Ordered trees,
             Cc0022: Finite Cartan types,
             Cc0023: Parking functions]

        """
        if isinstance(entry, FindStatCollection):
            return entry

        if isinstance(entry, (str, unicode)):
            # find by name in _findstat_collections
            for (id, c) in self._findstat_collections.iteritems():
                if entry.upper() in (c[0].upper(), c[1].upper(), c[2].upper()):
                    return self.element_class(self, id, c, None)

        elif isinstance(entry, (int, Integer)):
            # find by id in _findstat_collections
            for (id, c) in self._findstat_collections.iteritems():
                if entry == id:
                    return self.element_class(self, id, c, None)

        else:
            # find collection given an object or a constructor

            # unfortunately, we cannot test with
            # isinstance(_, SageObject), since this is True for
            # CartanType.

            # TODO: entry == c[4] will work rarely because c[4] might be a function!
            # also, the error handling is only necessary because of this...
            for (id, c) in self._findstat_collections.iteritems():
                try:
                    if isinstance(entry, c[3]) or entry == c[4]:
                        return self.element_class(self, id, c, None)

                except TypeError:
                    # examples are
                    # graphs:
                    # TypeError: cannot compare graph to non-graph (<class 'sage.combinat.permutation.Permutations'>)
                    # perfect matchings:
                    # TypeError: descriptor 'parent' of 'sage.structure.sage_object.SageObject' object needs an argument
                    pass

            # check whether entry is iterable (it's not a string!)
            try:
                obj = next(iter(entry))
                for (id, c) in self._findstat_collections.iteritems():
                    if isinstance(obj, c[3]):
                        return self.element_class(self, id, c, entry)

            except TypeError:
                pass

        raise ValueError("Could not find FindStat collection for %s." %str(entry))

    def _repr_(self):
        """
        Return the representation of the set of FindStat collections.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollections
            sage: FindStatCollections()                                         # optional -- internet
            Set of combinatorial collections used by FindStat
        """
        return "Set of combinatorial collections used by FindStat"

    def __iter__(self):
        """
        Return an iterator over all FindStat collections.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollections
            sage: [m for m in FindStatCollections()][0]                         # optional -- internet
            Cc0001: Permutations
        """
        for c in self._findstat_collections:
            yield FindStatCollection(c)

    Element = FindStatCollection

class FindStatMap(Element):
    r"""
    A FindStat map.

    :class:`FindStatMap` is a class representing a combinatorial
    map available in the FindStat database.

    The result of a :class:`findstat<FindStat>` query contains a
    (possibly empty) list of such maps.  This class provides methods
    to inspect various properties of these maps, in particular
    :meth:`code`.

    INPUT:

    - a string containing the FindStat name of the map, or an integer
      representing its FindStat id.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatMap
        sage: FindStatMap(71)                                                   # optional -- internet
        Mp00071: descent composition
        sage: FindStatMap("descent composition")                                # optional -- internet
        Mp00071: descent composition

    SEEALSO:

    :class:`FindStatMaps`

    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, entry):
        """
        Retrieve a map from the database.

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap("abcdefgh")                                       # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: Could not find FindStat map for abcdefgh.
        """
        return FindStatMaps()(entry)

    def __init__(self, parent, entry):
        """Initialize the map.

        This should only be called in
        :meth:`FindStatMaps()._element_constructor_` via
        `element_class`.

        INPUT:

        - ``parent`` -- :class:`FindStatMaps`.

        - ``entry`` -- a dictionary containing the properties of the
          map, such as its name, code, and so on.

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(62).parent()                                      # optional -- internet
            Set of combinatorial maps used by FindStat

        """
        self._map = entry
        Element.__init__(self, parent)

    def __reduce__(self):
        """Return a function and its arguments needed to create this
        map.

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: c = FindStatMap(62)                                           # optional -- internet
            sage: loads(dumps(c)) == c                                          # optional -- internet
            True

        """
        return (FindStatMap, (self.id(),))

    def id(self):
        r"""
        Return the FindStat identifier of the map.

        OUTPUT:

        The FindStat identifier of the map as an integer.

        EXAMPLES::

            sage: m = findstat("Permutations", lambda pi: pi.length())[1][1][0] # optional -- internet
            sage: m.id()                                                        # optional -- internet
            62
        """
        return self._map[FINDSTAT_MAP_IDENTIFIER]

    def id_str(self):
        r"""
        Return the FindStat identifier of the map.

        OUTPUT:

        The FindStat identifier of the map as a string.

        EXAMPLES::

            sage: m = findstat("Permutations", lambda pi: pi.length())[1][1][0] # optional -- internet
            sage: m.id_str()                                                    # optional -- internet
            'Mp00062'
        """

        id = str(self.id())
        return 'Mp00000'[:-len(id)] + id

    def _repr_(self):
        r"""
        Return the representation of the FindStat map.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(71)                                               # optional -- internet
            Mp00071: descent composition
        """
        return "%s: %s" %(self.id_str(), self._map[FINDSTAT_MAP_NAME])

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatMap, FindStatMaps
            sage: FindStatMap(71) == FindStatMap(71)                            # optional -- internet
            True

            sage: FindStatMap(62) == FindStatMap(71)                            # optional -- internet
            False

            sage: FindStatMap(71) != FindStatMap(71)                            # optional -- internet
            False

            sage: FindStatMap(62) != FindStatMap(71)                            # optional -- internet
            True

            sage: FindStatMap(62) == 1                                          # optional -- internet
            False

            sage: FindStatMap(62) != 1                                          # optional -- internet
            True

            sage: sorted(c for c in FindStatMaps())[0]                          # optional -- internet
            Mp00001: to semistandard tableau
        """
        return self.id().__cmp__(other.id())

    def name(self):
        r"""
        Return the FindStat name of the map.

        OUTPUT:

        The name of the map as a string, as used by FindStat.

        EXAMPLES::

            sage: m = findstat("Permutations", lambda pi: pi.length())[1][1][0] # optional -- internet
            sage: m.name()                                                      # optional -- internet
            u'inversion-number to major-index bijection'
        """
        return self._map[FINDSTAT_MAP_NAME]

    def description(self):
        r"""
        Return the FindStat description of the map.

        OUTPUT:

        The description as a string.

        EXAMPLES::

            sage: m = findstat("Permutations", lambda pi: pi.length())[1][1][0] # optional -- internet
            sage: print(m.description())                                    # optional -- internet,random
            Let $\sigma \in \mathcal{S}_n$ be a permutation.
            <BLANKLINE>
            Maps $\sigma$ to the permutation $\tau$ such that the major code of $\tau$ is given by the Lehmer code of $\sigma$.
            <BLANKLINE>
            In particular, the number of inversions of $\sigma$ equals the major index of $\tau$.
            <BLANKLINE>
            EXAMPLES:
            <BLANKLINE>
            $[3,4,1,2] \mapsto [3,1,4,2]$
        """
        return self._map[FINDSTAT_MAP_DESCRIPTION]

    def domain(self):
        r"""
        Return the FindStat collection which is the domain of the map.

        OUTPUT:

        The domain of the map as a :class:`FindStatCollection`.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: FindStatMap(71).domain()                                      # optional -- internet
            Cc0001: Permutations
        """
        return FindStatCollection(self._map[FINDSTAT_MAP_DOMAIN])

    def codomain(self):
        r"""
        Return the FindStat collection which is the codomain of the map.

        OUTPUT:

        The codomain of the map as a :class:`FindStatCollection`.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: FindStatMap(71).codomain()                                    # optional -- internet
            Cc0006: Integer compositions
        """
        return FindStatCollection(self._map[FINDSTAT_MAP_CODOMAIN])

    def code(self):
        r"""
        Return the code associated with the map.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: print(FindStatMap(71).code())                              # optional -- internet
            def descents_composition(elt):
                if len(elt) == 0:
                    return Composition([])
                d = [-1] + elt.descents() + [len(elt)-1]
                return Composition([ d[i+1]-d[i] for i in range(len(d)-1)])
        """
        return self._map[FINDSTAT_MAP_CODE]

    def code_name(self):
        r"""
        Return the name of the function defined by :meth:`code`.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: print(FindStatMap(71).code_name())                      # optional -- internet
            descents_composition
        """
        return self._map[FINDSTAT_MAP_CODE_NAME]

class FindStatMaps(Parent, UniqueRepresentation):
    r"""
    The class of FindStat maps.

    The elements of this class are combinatorial maps currently in
    FindStat.

    EXAMPLES:

    We can print a nice list of maps currently in FindStat, sorted by
    domain and codomain::

        sage: from sage.databases.findstat import FindStatMap, FindStatMaps
        sage: for m in sorted(FindStatMaps(), key=lambda m: (m.domain(), m.codomain)):    # optional -- internet,random
        ....:     print(m.domain().name().ljust(30)+" "+m.codomain().name().ljust(30)+" "+m.name())
        Permutation                    Standard tableau               Robinson-Schensted insertion tableau
        Permutation                    Integer partition              Robinson-Schensted tableau shape
        Permutation                    Binary tree                    to increasing tree
        ...

    """
    def __init__(self):
        """
        Fetch all the maps from FindStat.

        TESTS::

            sage: from sage.databases.findstat import FindStatMaps
            sage: M = FindStatMaps()                                            # optional -- internet
            sage: TestSuite(M).run()                                            # optional -- internet
        """
        self._findstat_maps = json.load(urlopen(FINDSTAT_URL_DOWNLOADS_MAPS))
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, entry):
        """Initialize a FindStat map.

        INPUT:

        - ``entry`` -- a string containing the FindStat name of the
          map, or an integer giving its id, or a dict containing all
          the information.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(71)                                               # optional -- internet
            Mp00071: descent composition

        """
        if isinstance(entry, FindStatMap):
            return entry

        elif entry in self._findstat_maps:
            return self.element_class(self, entry)

        elif isinstance(entry, (str, unicode)):
            # find by name in _findstat_maps
            for c in self._findstat_maps:
                if entry.upper() == c[FINDSTAT_MAP_NAME].upper():
                    return self.element_class(self, c)

        elif isinstance(entry, (int, Integer)):
            # find by id in _findstat_maps
            for c in self._findstat_maps:
                if entry == c[FINDSTAT_MAP_IDENTIFIER]:
                    return self.element_class(self, c)

        raise ValueError("Could not find FindStat map for %s." %str(entry))

    def _repr_(self):
        """
        Return the representation of the set of FindStat maps.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: FindStatMaps()                                                # optional -- internet
            Set of combinatorial maps used by FindStat
        """
        return "Set of combinatorial maps used by FindStat"

    def __iter__(self):
        """
        Return an iterator over all FindStat maps.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: [m for m in FindStatMaps()][0]                                # optional -- internet
            Mp00072: binary search tree: left to right

        """
        for m in self._findstat_maps:
            yield FindStatMap(m)

    Element = FindStatMap

findstat = FindStat()
