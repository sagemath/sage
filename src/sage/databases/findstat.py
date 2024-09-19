r"""
FindStat - the search engine for combinatorial statistics and maps

The interface to the FindStat database is ::

    sage: findstat()
    The Combinatorial Statistic Finder (https://www.findstat.org/)

We use the following three notions

- A *combinatorial collection* is a set `S` with interesting combinatorial properties,
- a *combinatorial map* is a combinatorially interesting map `f: S \to S'` between combinatorial collections, and
- a *combinatorial statistic* is a combinatorially interesting map `s: S \to \ZZ`.

You can use the FindStat interface to

- identify a combinatorial statistic or map given the values on a few small objects,
- obtain more terms, formulae, references, etc. for a given statistic or map,
- edit statistics and maps and submit new statistics.

The main entry points to the database are

- :func:`findstat` to search for matching statistics,
- :func:`findmap` to search for matching maps.

A guided tour
-------------

Retrieving information
^^^^^^^^^^^^^^^^^^^^^^

The most straightforward application of the FindStat interface is to
gather information about a combinatorial statistic.  To do this, we
supply :func:`findstat` with a list of ``(object, value)`` pairs.
For example::

    sage: PM = PerfectMatchings
    sage: r = findstat([(m, m.number_of_nestings()) for n in range(6) for m in PM(2*n)], depth=1); r # optional -- internet
    0: St000042oMp00116 (quality [100, 100])
    1: St000041 (quality [20, 100])
    ...

The result of this query is a list (presented as a
:class:`sage.databases.oeis.FancyTuple`) of matches.  Each match
consists of a :class:`FindStatCompoundStatistic` `s: S \to \ZZ` and
an indication of the quality of the match.

The precise meaning of the result is as follows:

    The composition `f_n \circ ... \circ f_2 \circ f_1` applied to
    the objects sent to FindStat agrees with all ``(object, value)``
    pairs of `s` in the database.  The optional parameter ``depth=1``
    limits the output to `n=1`.

    Suppose that the quality of the match is `(q_a, q_d)`.  Then
    `q_a` is the percentage of ``(object, value)`` pairs that are in
    the database among those which were sent to FindStat, and `q_d`
    is the percentage of ``(object, value)`` pairs with distinct
    values in the database among those which were sent to FindStat.

Put differently, if ``quality`` is not too small it is likely that
the statistic sent to FindStat equals `s \circ f_n \circ ... \circ
f_2 \circ f_1`.  If `q_a` is large, but `q_b` is small, then there
were many matches, but while the sought for statistic attains many
distinct values, the match found by FindStat covers only ``(object,
value)`` pairs for few values.

In the case at hand, for the match ``St000041``, the list of maps is
empty.  We can retrieve the description of the statistic from the
database as follows::

    sage: print(r[1].statistic().description())                                 # optional -- internet
    The number of nestings of a perfect matching.
    <BLANKLINE>
    <BLANKLINE>
    This is the number of pairs of edges $((a,b), (c,d))$ such that $a\le c\le d\le b$. i.e., the edge $(c,d)$ is nested inside $(a,b)$...

We can check the references::

    sage: r[1].statistic().references()                                         # optional -- internet
    0: [1]  de Médicis, A., Viennot, X. G., Moments des $q$-polynômes de Laguerre et la bijection de Foata-Zeilberger [[MathSciNet:1288802]]
    1: [2]  Simion, R., Stanton, D., Octabasic Laguerre polynomials and permutation statistics [[MathSciNet:1418763]]...

If you prefer, you can look at this information also in your browser::

    sage: r[1].statistic().browse()                                             # optional -- webbrowser

Another interesting possibility is to look for equidistributed
statistics.  Instead of submitting a list of ``(object, value)``
pairs, we pass a list of pairs ``(objects, values)``::

    sage: data = [(PM(2*n), [m.number_of_nestings() for m in PM(2*n)]) for n in range(5)]
    sage: findstat(data, depth=0)                                               # optional -- internet
    0: St000041 (quality [99, 100])
    1: St000042 (quality [99, 100])

This results tells us that the database contains another entry that
is equidistributed with the number of nestings on perfect matchings
of size at most `10`, namely the number of crossings.  Note that
there is a limit on the number of elements FindStat accepts for a
query, which is currently `1000`.  Queries with more than `1000`
elements are truncated.

Let us now look at a slightly more complicated example, where the
submitted statistic is the composition of a sequence of combinatorial
maps and a statistic known to FindStat.  We use the occasion to
advertise yet another way to pass values to FindStat::

    sage: r = findstat(Permutations, lambda pi: pi.saliances()[0], depth=2); r  # optional -- internet
    0: St000740oMp00066 with offset 1 (quality [100, 100])
    1: St000740oMp00087 with offset 1 (quality [100, 100])
    ...

Note that some of the matches are up to a global offset.  For
example, we have::

    sage: r[1].info()                                                           # optional -- internet
    after adding 1 to every value
    and applying
        Mp00087: inverse first fundamental transformation: Permutations -> Permutations
    to the objects (see `.compound_map()` for details)
    <BLANKLINE>
    your input matches
        St000740: The last entry of a permutation.
    <BLANKLINE>
    among the values you sent, 100 percent are actually in the database,
    among the distinct values you sent, 100 percent are actually in the database

Let us pick another particular result::

    sage: s = findstat("St000051oMp00061oMp00069")                              # optional -- internet
    sage: s.info()                                                              # optional -- internet
        Mp00069: complement: Permutations -> Permutations
        Mp00061: to increasing tree: Permutations -> Binary trees
        St000051: The size of the left subtree of a binary tree.

To obtain the value of the statistic sent to FindStat on a given
object, apply the maps in the list in the given order to this object,
and evaluate the statistic on the result.  For example, let us check
that the result given by FindStat agrees with our statistic on the
following permutation::

    sage: pi = Permutation([3,1,4,5,2]); pi.saliances()[0]
    3

We first have to find out, what the maps and the statistic actually do::

    sage: print(s.statistic().description())                                    # optional -- internet
    The size of the left subtree of a binary tree.

    sage: print(s.statistic().sage_code())                                      # optional -- internet
    def statistic(T):
        return T[0].node_number()

    sage: print("\n\n".join(m.sage_code() for m in s.compound_map()))           # optional -- internet
    def mapping(sigma):
        return sigma.complement()
    <BLANKLINE>
    def mapping(sigma):
        return sigma.increasing_tree_shape()

So, the following should coincide with what we sent FindStat::

    sage: pi.complement().increasing_tree_shape()[0].node_number()
    3

Editing and submitting statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Of course, often a statistic will not be in the database::

    sage: s = findstat([(d, randint(1,1000)) for d in DyckWords(4)]); s         # optional -- internet
    St000000: a new statistic on Dyck paths

In this case, and if the statistic might be "interesting", please
consider submitting it to the database using
:meth:`FindStatStatistic.submit`.

Also, you may notice omissions, typos or even mistakes in the
description, the code and the references.  In this case, simply
replace the value by using :meth:`FindStatFunction.set_description`,
:meth:`FindStatStatistic.set_code` or
:meth:`FindStatFunction.set_references_raw`, and then
:meth:`FindStatStatistic.submit` your changes for review by the
FindStat team.

AUTHORS:

- Martin Rubey (2015): initial version
- Martin Rubey (2020): rewrite, adapt to new FindStat API
"""
# ****************************************************************************
#       Copyright (C) 2015 Martin Rubey <martin.rubey@tuwien.ac.at>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.lazy_list import lazy_list
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.sets_cat import Sets
from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp

from sage.misc.verbose import verbose
from sage.rings.integer import Integer
from sage.databases.oeis import FancyTuple

from ast import literal_eval
from copy import deepcopy
import re
import webbrowser
import tempfile
import inspect
import html
import requests
from json.decoder import JSONDecodeError
import itertools

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
from sage.combinat.parking_functions import ParkingFunction, ParkingFunctions
from sage.combinat.perfect_matching import PerfectMatching, PerfectMatchings
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.posets.posets import Poset, FinitePoset
from sage.combinat.posets.lattices import LatticePoset, FiniteLatticePoset
from sage.combinat.posets.poset_examples import Posets
from sage.combinat.tableau import SemistandardTableau, SemistandardTableaux, StandardTableau, StandardTableaux
from sage.combinat.set_partition import SetPartition, SetPartitions
from sage.combinat.skew_partition import SkewPartition, SkewPartitions
from sage.graphs.graph_generators import graphs
from sage.combinat.words.word import Word
from sage.combinat.words.words import Words
from sage.combinat.words.abstract_word import Word_class
from sage.combinat.colored_permutations import SignedPermutation, SignedPermutations
from sage.combinat.plane_partition import PlanePartition
from sage.combinat.decorated_permutation import DecoratedPermutation, DecoratedPermutations
from sage.combinat.set_partition_ordered import OrderedSetPartition, OrderedSetPartitions

######################################################################
# the FindStat URLs
FINDSTAT_URL = 'https://www.findstat.org/'
FINDSTAT_API = FINDSTAT_URL + "api/"
FINDSTAT_API_COLLECTIONS = FINDSTAT_API + 'CollectionsDatabase/'
FINDSTAT_API_STATISTICS = FINDSTAT_API + 'StatisticsDatabase/'
FINDSTAT_API_MAPS = FINDSTAT_API + 'MapsDatabase/'

FINDSTAT_URL_LOGIN = FINDSTAT_URL + "?action=login"
FINDSTAT_URL_COLLECTIONS = FINDSTAT_URL + 'CollectionsDatabase/'
FINDSTAT_URL_STATISTICS = FINDSTAT_URL + 'StatisticsDatabase/'
FINDSTAT_URL_EDIT_STATISTIC = FINDSTAT_URL + 'EditStatistic/'
FINDSTAT_URL_NEW_STATISTIC = FINDSTAT_URL + 'NewStatistic/'
FINDSTAT_URL_MAPS = FINDSTAT_URL + 'MapsDatabase/'
FINDSTAT_URL_EDIT_MAP = FINDSTAT_URL + 'EditMap/'
FINDSTAT_URL_NEW_MAP = FINDSTAT_URL + 'NewMap/'

######################################################################
# the number of values FindStat allows to search for at most
FINDSTAT_MAX_VALUES = 1000
# the number of values FindStat needs at least to search for
FINDSTAT_MIN_VALUES = 3
# the maximal number of maps that FindStat composes by default to find a match
FINDSTAT_DEFAULT_DEPTH = 2
# the number of values FindStat allows to submit at most
FINDSTAT_MAX_SUBMISSION_VALUES = 1200

# separates name from description
FINDSTAT_SEPARATOR_NAME = "\n"
# separates values
FINDSTAT_VALUE_SEPARATOR = ";"
FINDSTAT_MAP_SEPARATOR = "o"
# regexp to recognize a statistic identifier
FINDSTAT_STATISTIC_REGEXP = '^St[0-9]{6}$'
FINDSTAT_MAP_REGEXP = '^Mp[0-9]{5}$'
FINDSTAT_COLLECTION_REGEXP = '^Cc[0-9]{4}$'
FINDSTAT_STATISTIC_PADDED_IDENTIFIER = "St%06d"
FINDSTAT_MAP_PADDED_IDENTIFIER = "Mp%05d"
FINDSTAT_COLLECTION_PADDED_IDENTIFIER = "Cc%04d"

######################################################################

# the format string for using POST
# WARNING: we use html.escape to avoid injection problems, thus we expect double quotes as field delimiters.
FINDSTAT_POST_HEADER = """
<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.3.2/jquery.min.js"></script>
<script>
    $(document).ready(function() {$("#form").submit(); });
</script>
"""
FINDSTAT_NEWSTATISTIC_FORM_HEADER = '<form id="form" name="NewStatistic" action="%s" enctype="multipart/form-data" method="post" />'
FINDSTAT_NEWMAP_FORM_HEADER = '<form id="form" name="NewMap" action="%s" enctype="multipart/form-data" method="post" />'
FINDSTAT_FORM_FORMAT = '<input type="hidden" name="%s" value="%s" />'
FINDSTAT_FORM_FOOTER = '</form>'

######################################################################


class FindStat(UniqueRepresentation, SageObject):
    r"""
    The Combinatorial Statistic Finder.

    :class:`FindStat` is a class preserving user information.
    """
    def __init__(self):
        r"""
        Initialize the database.

        TESTS::

            sage: findstat()
            The Combinatorial Statistic Finder (https://www.findstat.org/)
        """
        # user credentials if provided
        self._user_name = ""
        self._user_email = ""
        self._allow_execution = False

    def __repr__(self):
        r"""
        Return the representation of ``self``.

        EXAMPLES::

            sage: findstat()
            The Combinatorial Statistic Finder (https://www.findstat.org/)
        """
        return "The Combinatorial Statistic Finder (%s)" % FINDSTAT_URL

    def browse(self):
        r"""
        Open the FindStat web page in a browser.

        EXAMPLES::

            sage: findstat().browse()                                           # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL)

    def set_user(self, name=None, email=None):
        r"""
        Set the user for the session.

        INPUT:

        - ``name`` -- the name of the user

        - ``email`` -- an email address of the user

        This information is used when submitting a statistic with
        :meth:`FindStatStatistic.submit`.

        EXAMPLES::

            sage: findstat().set_user(name='Anonymous', email='invalid@org')

        .. NOTE::

            It is usually more convenient to login into the FindStat
            web page using the :meth:`login` method.
        """
        if not isinstance(name, str):
            raise ValueError("the given name (%s) should be a string" % name)
        if not isinstance(email, str):
            raise ValueError("the given email (%s) should be a string" % email)
        self._user_name = name
        self._user_email = email

    def user_name(self):
        """
        Return the user name used for submissions.

        EXAMPLES::

            sage: findstat().set_user(name='Anonymous', email='invalid@org')
            sage: findstat().user_name()
            'Anonymous'
        """
        return self._user_name

    def user_email(self):
        """
        Return the user name used for submissions.

        EXAMPLES::

            sage: findstat().set_user(name='Anonymous', email='invalid@org')
            sage: findstat().user_email()
            'invalid@org'
        """
        return self._user_email

    def login(self):
        r"""
        Open the FindStat login page in a browser.

        EXAMPLES::

            sage: findstat().login()                                            # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_LOGIN)


######################################################################
# tools
######################################################################
def _get_json(url, **kwargs):
    """
    Return the json response or raise an error.

    EXAMPLES::

        sage: from sage.databases.findstat import _get_json, FINDSTAT_API_MAPS
        sage: _get_json(FINDSTAT_API_MAPS + "?fields=yyy")                      # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: E018: Unknown fields in Map, Combinatorial map: to semistandard tableau via monotone triangles: ['yyy']
    """
    response = requests.get(url)
    if response.ok:
        try:
            result = response.json(**kwargs)
        except JSONDecodeError:
            raise ValueError(response.text)
        if "error" in result:
            raise ValueError(result["error"])
        return result
    raise ConnectionError(response.text)


def _post_json(url, data, **kwargs):
    """
    Return the json response or raise an error.

    EXAMPLES::

        sage: from sage.databases.findstat import _post_json, FINDSTAT_API_STATISTICS
        sage: _post_json(FINDSTAT_API_STATISTICS, {"fields": "yyy"})            # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: E018: Unknown fields in Statistic, Combinatorial statistic: St000001: ['yyy']
    """
    response = requests.post(url, data=data)
    if response.ok:
        try:
            result = response.json(**kwargs)
        except JSONDecodeError:
            raise ValueError(response.text)
        if "error" in result:
            raise ValueError(result["error"])
        return result
    raise ConnectionError(response.text)


def _submit(args, url):
    r"""
    Open a post form containing fields for each of the arguments,
    which is sent to the given url.

    INPUT:

    - args, a dictionary whose keys are the form fields

    - url, the goal of the post request

    EXAMPLES::

        sage: from sage.databases.findstat import _submit, FINDSTAT_NEWSTATISTIC_FORM_HEADER, FINDSTAT_URL_NEW_STATISTIC
        sage: url = FINDSTAT_NEWSTATISTIC_FORM_HEADER % FINDSTAT_URL_NEW_STATISTIC
        sage: args = {"OriginalStatistic": "St000000",
        ....:         "Domain": "Cc0001",
        ....:         "Values": "[1,3,2]=>17\n[1,2,3]=>83",
        ....:         "Description": "Not a good statistic",
        ....:         "References": "[[arXiv:1102.4226]]",
        ....:         "Code": "def statistic(p):\n return 1",
        ....:         "SageCode": "def statistic(p):\n return 1",
        ....:         "CurrentAuthor": "",
        ....:         "CurrentEmail": ""}
        sage: _submit(args, url)                                                # optional -- webbrowser
    """
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.html', delete=False)
    verbose("Created temporary file %s" % f.name, caller_name='FindStat')
    f.write(FINDSTAT_POST_HEADER)
    f.write(url)
    for key, value in args.items():
        if value:
            verbose("writing argument %s" % key, caller_name='FindStat')
            value_encoded = html.escape(value, quote=True)
            html_content = FINDSTAT_FORM_FORMAT % (key, value_encoded)
            f.write(html_content)
        else:
            verbose("skipping argument %s because it is empty" % key, caller_name='FindStat')
    f.write(FINDSTAT_FORM_FOOTER)
    f.close()
    verbose("Opening file with webbrowser", caller_name='FindStat')
    webbrowser.open(f.name)


def _data_to_str(data, domain, codomain=None):
    """
    Return a string representation of the given list of ``(objects,
    values)`` pairs suitable for a FindStat query.

    INPUT:

    - ``data`` -- list of lists of objects

    - ``domain`` -- :class:`FindStatCollection`

    - ``codomain`` -- (optional) :class:`FindStatCollection` or ``None``

    If ``codomain`` is ``None``, the values are treated as integers.

    TESTS::

        sage: n = 3; l = lambda i: [pi for pi in Permutations(n) if pi(1) == i]
        sage: data = [([pi for pi in l(i)], [pi(1) for pi in l(i)]) for i in range(1,n+1)]
        sage: data.append(([Permutation([1,2])], [1]))
        sage: from sage.databases.findstat import FindStatCollection, _data_to_str
        sage: print(_data_to_str(data, FindStatCollection(1)))                  # optional -- internet
        [1, 2, 3]
        [1, 3, 2]
        ====> 1;1
        [2, 1, 3]
        [2, 3, 1]
        ====> 2;2
        [3, 1, 2]
        [3, 2, 1]
        ====> 3;3
        [1, 2]
        ====> 1
    """
    to_str_dom = domain.to_string()
    if codomain is None:
        to_str_codom = str
    else:
        to_str_codom = codomain.to_string()

    return "\n".join("\n".join(to_str_dom(element) for element in elements)
                     + "\n====> "
                     + FINDSTAT_VALUE_SEPARATOR.join(to_str_codom(value)
                                                     for value in values)
                     for elements, values in data)


def _data_from_iterable(iterable, mapping=False, domain=None,
                        codomain=None, check=True):
    """
    Return a list of pairs of lists of the same size, domain, and if
    applicable, codomain.

    INPUT:

    - iterable, a pair of lists of the same size, or an iterable of
      pairs, such that every pair consists of two iterables of the
      same size, or a single element and a single value.  Every
      object must be a :class:`SageObject`.

    - ``mapping`` -- boolean (default: ``False``); ``False``, if the codomain
      is ``Integer`` and ``True`` if it is a FindStat collection

    - ``domain`` -- (default: ``None``) a :class:`FindStatCollection`, if
      ``None`` it is guessed from the iterable

    - ``codomain`` -- (default: ``None``) a :class:`FindStatCollection`, if
      ``None`` it is guessed from the iterable

    TESTS::

        sage: from sage.databases.findstat import FindStatCollection, _data_from_iterable
        sage: PM = PerfectMatchings
        sage: l = [(PM(2*n), [m.number_of_nestings() for m in PM(2*n)]) for n in range(5)]
        sage: _data_from_iterable(l)                                            # optional -- internet
        (lazy list [([[]], [0]), ([[(1, 2)]], [0]), ([[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]], [0, 0, 1]), ...],
         a subset of Cc0012: Perfect matchings)
        sage: domain = FindStatCollection("Set Partitions")                     # optional -- internet
        sage: _data_from_iterable(l, domain=domain)                             # optional -- internet
        (lazy list [([[]], [0]), ([[(1, 2)]], [0]), ([[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]], [0, 0, 1]), ...],
         Cc0009: Set partitions)

    Check that passing a list with a single ``(elements, values)`` pair works::

        sage: l = [(PM(4), [0, 0, 1])]
        sage: _data_from_iterable(l, domain=domain)                             # optional -- internet
        (lazy list [([[(1, 2), (3, 4)], [(1, 3), (2, 4)], [(1, 4), (2, 3)]], [0, 0, 1])],
         Cc0009: Set partitions)
    """
    if isinstance(iterable, dict):
        iterator = iter(iterable.items())
    else:
        iterator = iter(iterable)
    query0 = next(iterator)
    try:
        query1 = next(iterator)
    except StopIteration:
        # query0 == (elts, vals), which we interpret as [(elts, vals)]
        query0, query1 = query0[0], query0[1]
        # next(iterator) will raise StopIteration
    try:
        query2 = next(iterator)
        pre_data = [query0, query1, query2]
    except StopIteration:
        # (query0, query1) == (elts, vals) or
        # (query0, query1) == [(elts, vals), (elts, vals)]
        # in the former case, elts has at least 3 elements
        # in both cases, query0 and query1 are not objects
        # if query0 is a long list, we have to raise an error anyway
        elts, vals = list(query0), list(query1)
        if len(elts) == 2:
            if len(vals) != 2:
                raise ValueError("cannot interpret the given argument as a FindStat query")
            pre_data = [elts, vals]
        else:
            if len(elts) != len(vals):
                raise ValueError("FindStat expects the same number of objects (got %s) as values (got %s)" % (len(elts), len(vals)))
            pre_data = [(elts, vals)]

    # pre_data is a list of all elements of the iterator accessed so
    # far, for each of its elements and also the remainder ot the
    # iterator, each element is either a pair ``(object, value)`` or
    # a pair ``(objects, values)``
    elts, vals = pre_data[0]
    if domain is None:
        domain = FindStatCollection(elts)
    if mapping and codomain is None:
        codomain = FindStatCollection(vals)

    all_elements = set()

    def sanitize_pair(elts, vals):
        if domain.is_element(elts):
            elts = [elts]
            if mapping:
                vals = [vals]
            else:
                vals = [Integer(vals)]
        else:
            elts = list(elts)
            if check:
                bad = [elt for elt in elts if not domain.is_element(elt)]
                assert not bad, "%s are not elements of %s" % (bad, domain)
            if mapping:
                vals = list(vals)
            else:
                vals = list(map(Integer, vals))
        if len(elts) != len(vals):
            raise ValueError("FindStat expects the same number of objects as values in each pair")
        if check and mapping:
            bad = [elt for elt in vals if not codomain.is_element(elt)]
            assert not bad, "%s are not elements of %s" % (bad, codomain)
        for elt in elts:
            if elt in all_elements:
                raise ValueError("FindStat expects that every object occurs at most once: %s" % elt)
            all_elements.add(elt)

        return elts, vals

    lazy_data = lazy_list((sanitize_pair(elts, vals)
                           for elts, vals in iterator),
                          initial_values=[sanitize_pair(elts, vals)
                                          for elts, vals in pre_data])
    if mapping:
        return lazy_data, domain, codomain
    return lazy_data, domain


def _data_from_function(function, domain):
    """
    Return a lazy list of pairs of singleton lists of the same size,
    domain.

    INPUT:

    - ``function`` -- a callable
    - ``domain`` -- a :class:`FindStatCollection`

    If ``function`` returns the value ``None``, the pair is omitted.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatCollection, _data_from_function
        sage: domain = FindStatCollection(1)                                    # optional -- internet
        sage: _data_from_function(lambda pi: pi[0], domain)                     # optional -- internet
        lazy list [([[1]], [1]), ([[1, 2]], [1]), ([[2, 1]], [2]), ...]
    """
    return lazy_list(([elt], [value])
                     for elt, value in domain.first_terms(function)
                     if value is not None)


def _data_from_data(data, max_values):
    """
    Return the first few pairs (of lists of the same size) with a
    total of at most ``max_values`` objects in the range of the
    collection.

    INPUT:

    - ``data`` -- an iterable over pairs of lists of the same size

    - ``max_values`` -- the maximal number of objects (and values) to
      return

    We assume that the number of elements in each pair weakly
    increases, to decide when to stop.

    TESTS::

        sage: from sage.databases.findstat import _data_from_data
        sage: PM = PerfectMatchings
        sage: l = [(PM(2*n), (m.number_of_nestings() for m in PM(2*n))) for n in range(2,100)]
        sage: [(a, list(b)) for a, b in _data_from_data(l, 35)]
        [(Perfect matchings of {1, 2, 3, 4}, [0, 0, 1]),
         (Perfect matchings of {1, 2, 3, 4, 5, 6},
          [0, 0, 1, 1, 2, 2, 1, 0, 0, 0, 1, 1, 1, 2, 3])]
    """
    query = []
    total = min(max_values, FINDSTAT_MAX_VALUES)
    iterator = iter(data)
    while total > 0:
        try:
            elts, vals = next(iterator)
        except StopIteration:
            break
        if total >= len(elts):
            query.append((elts, vals))
            total -= len(elts)
        else:
            break # assuming that the next pair is even larger

    return query


def _distribution_from_data(data, domain, max_values, generating_functions=False):
    """
    Return the first few pairs (of lists of the same size) with a
    total of at most ``max_values`` objects in the range of the
    collection, combined by level.

    INPUT:

    - ``data`` -- an iterable over pairs of lists of the same size

    - ``domain`` -- a :class:`FindStatCollection`

    - ``max_values`` -- the maximal number of objects (and values) to
      return

    TESTS::

        sage: from sage.databases.findstat import _distribution_from_data, FindStatCollection
        sage: n = 3; l = lambda i: [pi for pi in Permutations(n) if pi(1) == i]
        sage: data = [([pi for pi in l(i)], [pi(1) for pi in l(i)]) for i in range(1,n+1)]
        sage: cc = FindStatCollection(1)                                        # optional -- internet
        sage: _distribution_from_data(data, cc, 10)                             # optional -- internet
        [([[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]],
          [1, 1, 2, 2, 3, 3])]
    """
    lvl_dict = {} # lvl: elts, vals
    total = max_values
    iterator = iter(data)
    levels_with_sizes = domain.levels_with_sizes()
    while total > 0:
        try:
            elts, vals = next(iterator)
        except StopIteration:
            break
        if total < len(elts):
            break
        lvl = domain.element_level(elts[0])
        if generating_functions and lvl not in levels_with_sizes:
            continue
        if levels_with_sizes[lvl] > total:
            # we assume that from now on levels become even larger
            break
        if not all(domain.element_level(elt) == lvl for elt in elts[1:]):
            raise ValueError("cannot combine %s into a distribution" % elts)
        lvl_elts, lvl_vals = lvl_dict.get(lvl, [[], []])
        lvl_dict[lvl] = (lvl_elts + elts, lvl_vals + vals)
        if levels_with_sizes[lvl] == len(lvl_dict[lvl][0]):
            total -= levels_with_sizes[lvl]

    if generating_functions:
        return {lvl: {val: vals.count(val) for val in set(vals)}
                for lvl, (elts, vals) in lvl_dict.items()
                if levels_with_sizes[lvl] == len(vals)}

    return [(elts, vals) for lvl, (elts, vals) in lvl_dict.items()
            if levels_with_sizes[lvl] == len(elts)]


def _generating_functions_from_dict(gfs, style):
    """
    Convert the generating functions given as a dictionary into a
    desired style.

    INPUT:

    - ``gfs`` -- dictionary whose keys are the levels and whose values
      are dictionaries from values to multiplicities

    - ``style`` -- one of ``'dictionary'``, ``'list'`` or
      ``'polynomial'``

    .. SEEALSO::

        - :meth:`FindStatCombinatorialStatistic.generating_functions()`

    TESTS::

        sage: from sage.databases.findstat import _generating_functions_from_dict
        sage: d = {n: {2*i: n*i+1 for i in range(n, 2*n+1)} for n in range(3)}; d
        {0: {0: 1}, 1: {2: 2, 4: 3}, 2: {4: 5, 6: 7, 8: 9}}

        sage: _generating_functions_from_dict(d, "polynomial")
        {0: 1, 1: 3*q^4 + 2*q^2, 2: 9*q^8 + 7*q^6 + 5*q^4}

        sage: _generating_functions_from_dict(d, "list")
        {0: [1], 1: [2, 0, 3], 2: [5, 0, 7, 0, 9]}
    """
    if style == "dictionary":
        return gfs
    if style == "list":
        return {level: [gen_dict.get(deg, 0)
                        for deg in range(min(gen_dict),
                                         max(gen_dict)+1)]
                for level, gen_dict in gfs.items() if gen_dict}
    if style == "polynomial":
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.integer_ring import ZZ
        P = PolynomialRing(ZZ, "q", sparse=True)
        q = P.gen()
        return {level: sum(coefficient * q**exponent
                           for exponent, coefficient in gen_dict.items())
                for level, gen_dict in gfs.items()}

    raise ValueError("the argument 'style' (='%s') must be 'dictionary', 'polynomial', or 'list'" % style)


def _get_code_from_callable(function):
    """
    Return code given a callable, if possible.

    TESTS::

        sage: from sage.databases.findstat import _get_code_from_callable
        sage: @cached_function
        ....: def statistic(pi):
        ....:     pi = Permutations(len(pi))(pi)
        ....:     if pi.is_one():
        ....:         return 1
        ....:     return sum(statistic(pi.apply_simple_reflection_right(i))
        ....:                for i in pi.descents())
        sage: print(_get_code_from_callable(statistic))
        @cached_function
        def statistic(pi):
            pi = Permutations(len(pi))(pi)
            if pi.is_one():
                return Integer(1)
            return sum(statistic(pi.apply_simple_reflection_right(i))
                       for i in pi.descents())
    """
    code = ""
    if function is not None:
        from sage.misc.cachefunc import CachedFunction
        try:
            if isinstance(function, CachedFunction):
                code = inspect.getsource(function.f)
            else:
                code = inspect.getsource(function)
        except (OSError, TypeError):
            verbose("inspect.getsource could not get code from function provided",
                    caller_name='FindStat')
    return code


def findstat(query=None, values=None, distribution=None, domain=None,
             depth=FINDSTAT_DEFAULT_DEPTH, max_values=FINDSTAT_MAX_VALUES):
    r"""
    Return matching statistics.

    INPUT:

    One of the following:

    - an integer or a string representing a valid FindStat identifier
      (e.g. 45 or 'St000045').  The keyword arguments ``depth`` and
      ``max_values`` are ignored, ``values`` and ``distribution``
      must be ``None``.

    - a list of pairs of the form ``(object, value)``, or a
      dictionary from Sage objects to integer values.  The keyword
      arguments ``depth`` and ``max_values`` are passed to the
      finder, ``values`` and ``distribution`` must be ``None``.

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

    - a collection and a callable.  The callable is used to generate
      ``max_values`` ``(object, value)`` pairs.  The number of terms
      generated may also be controlled by passing an iterable
      collection, such as ``Permutations(3)``.  The keyword arguments
      ``depth`` and ``max_values`` are passed to the finder.

    OUTPUT: an instance of a :class:`FindStatStatistic`, represented by

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

        sage: findstat('St000041')                                              # optional -- internet
        St000041: The number of nestings of a perfect matching.

        sage: findstat(51)                                                      # optional -- internet
        St000051: The size of the left subtree of a binary tree.

        sage: findstat('St000042oMp00116')                                      # optional -- internet
        St000042oMp00116

    The database can be searched by providing a list of pairs::

        sage: l = [m for n in range(1, 4) for m in PerfectMatchings(2*n)]
        sage: findstat([(m, m.number_of_nestings()) for m in l], depth=0)       # optional -- internet
        0: St000041 (quality [100, 100])

    or a dictionary::

        sage: findstat({m: m.number_of_nestings() for m in l}, depth=0)         # optional -- internet
        0: St000041 (quality [100, 100])

    Note however, that the results of these two queries need not
    compare equal, because we compare queries by the data
    sent, and the ordering of the data might be different.

    Another possibility is to send a collection and a function.  In
    this case, the function is applied to the first few objects of
    the collection::

        sage: findstat("Perfect Matchings", lambda m: m.number_of_nestings(), depth=0)    # optional -- internet
        0: St000041 (quality [20, 100])

    To search for a distribution, send a list of lists, or a single pair::

        sage: PM = PerfectMatchings(10); findstat((PM, [m.number_of_nestings() for m in PM]), depth=0) # optional -- internet
        0: St000042 (quality [100, 100])
        1: St000041 (quality [9, 100])

    Alternatively, specify the ``distribution`` parameter::

        sage: findstat(12, distribution=lambda m: m.number_of_nestings(), depth=0) # optional -- internet
        0: St000041 (quality [100, 100])
        1: St000042 (quality [100, 100])

    Note that there is a limit, ``FINDSTAT_MAX_VALUES``, on the number
    of elements that may be submitted to FindStat, which is currently
    1000.  Therefore, the interface tries to truncate queries
    appropriately, but this may be impossible, especially with
    distribution searches::

        sage: PM = PerfectMatchings(12); PM.cardinality()                       # optional -- internet
        10395
        sage: findstat((PM, [1 for m in PM]))                                   # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: E016: The statistic finder was unable to perform a search on your data. The following errors have occured:
        <BLANKLINE>
        You passed too few elements (0 < 3) to FindStat!

    Finally, we can also retrieve all statistics with a given domain::

        sage: findstat("Cc0024")                                                # optional -- internet
        Set of combinatorial statistics with domain Cc0024: Binary words in FindStat

        sage: findstat(domain='Cores')                                          # optional -- internet
        Set of combinatorial statistics with domain Cc0013: Cores in FindStat

    TESTS::

        sage: findstat("Permutations", lambda x: 1, depth='x')                  # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: E021: Depth should be a nonnegative integer at most 9, but is x.

        sage: findstat("Permutations", lambda x: 1, depth=100)                  # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: E021: Depth should be a nonnegative integer at most 9, but is 100.

        sage: S = Permutation
        sage: findstat([(S([1,2]), 1), ([S([1,3,2]), S([1,2])], [2,3])])        # optional -- internet
        Traceback (most recent call last):
        ...
        ValueError: FindStat expects that every object occurs at most once: [1, 2]

    Check that values which can be converted to integers are supported::

        sage: findstat([(m, m.number_of_nestings()/1) for m in PerfectMatchings(10)], depth=0) # optional -- internet
        0: St000041 (quality [9, 100])

    Check that ``None`` values are omitted::

        sage: findstat("graphs", lambda g: g.diameter() if g.is_connected() else None, max_values=100, depth=0) # optional -- internet
        0: St000259 (quality [100, 100])
    """
    try:
        max_values = int(max_values)
        assert 0 <= max_values <= FINDSTAT_MAX_VALUES
    except (ValueError, AssertionError):
        raise ValueError("the maximal number of values for a FindStat query must be a nonnegative integer less than or equal to %i" % FINDSTAT_MAX_VALUES)

    check_collection = True

    def get_values(raw, domain=None):
        if callable(raw):
            known_terms = _data_from_function(raw, domain)
            function = raw
        else:
            known_terms, domain = _data_from_iterable(raw, domain=domain,
                                                      mapping=False,
                                                      check=check_collection)
            function = None
        data = _data_from_data(known_terms, max_values)
        return known_terms, data, domain, function

    def get_distribution(raw, domain=None):
        if callable(raw):
            known_terms = _data_from_function(raw, domain)
            function = raw
        else:
            known_terms, domain = _data_from_iterable(raw, domain=domain,
                                                      mapping=False,
                                                      check=check_collection)
            function = None
        data = _distribution_from_data(known_terms, domain, max_values)
        return known_terms, data, domain, function

    ######################################################################
    if query is None and values is None and distribution is None and domain is None:
        return FindStat()

    if values is not None and distribution is not None:
        raise ValueError("not both of `values` and `distribution` may be given for a FindStat query")

    if values is None and distribution is None:
        if query is None:
            return FindStatStatistics(domain=domain)

        if domain is None:
            if isinstance(query, (int, Integer, FindStatCombinatorialStatistic)):
                return FindStatStatistic(query)

            if isinstance(query, str):
                if re.match(FINDSTAT_COLLECTION_REGEXP, query):
                    return FindStatStatistics(domain=query)

                return FindStatStatistic(query)

        values, query = query, None

    if query is not None:
        if domain is None:
            domain, query = query, None
        else:
            raise ValueError("the given arguments cannot be used for a FindStat search")

    if domain is not None:
        domain = FindStatCollection(domain)

    if values is not None:
        if isinstance(values, (int, Integer, str, FindStatCombinatorialStatistic)):
            if domain is not None:
                raise ValueError("the domain must not be provided if a statistic identifier is given")
            return FindStatStatisticQuery(values_of=values, depth=depth)

        known_terms, data, domain, function = get_values(values, domain)
        return FindStatStatisticQuery(data=data, domain=domain, depth=depth,
                                      known_terms=known_terms, function=function)

    if distribution is not None:
        if isinstance(distribution, (int, Integer, str, FindStatCombinatorialStatistic)):
            if domain is not None:
                raise ValueError("the domain must not be provided if a statistic identifier is given")
            return FindStatStatisticQuery(distribution_of=distribution, depth=depth)

        known_terms, data, domain, function = get_distribution(distribution, domain)
        return FindStatStatisticQuery(data=data, domain=domain, depth=depth,
                                      known_terms=known_terms, function=function)

    raise ValueError("the given arguments cannot be used for a FindStat search")


def findmap(*args, **kwargs):
    r"""
    Return matching maps.

    INPUT:

    Valid keywords are: ``domain``, ``codomain``, ``values``,
    ``distribution``, ``depth`` and ``max_values``. They have the
    following meanings:

    - ``depth`` -- nonnegative integer (default: ``FINDSTAT_DEFAULT_DEPTH``);
      specifying how many maps to apply to generate the given map

    - ``max_values`` -- integer (default: ``FINDSTAT_MAX_VALUES``); specifying
      how many values are sent to the finder

    - ``domain``, ``codomain`` -- integer or string of the form
      ``Cc1234``, designates the domain and codomain of the sought
      for maps

    - ``values``, ``distribution`` -- data specifying the values or
      distribution of values of the sought for maps.  The keyword
      arguments ``depth`` and ``max_values`` are passed to the
      finder.  The data may be specified in one of the following
      forms:

        - a list of pairs of the form ``(object, value)``, or a
          dictionary from Sage objects to Sage objects.

        - a list of pairs of the form ``(list of objects, list of
          values)``, or a single pair of the form ``(list of objects,
          list of values)``.  In each pair there should be as many
          objects as values.

        - a callable.  In this case, the domain must be specified,
          also.  The callable is then used to generate ``max_values``
          ``(object, value)`` pairs.

          The number of terms generated may also be controlled by
          passing an iterable collection, such as
          ``Permutations(3)``.

    ``findmap`` also accepts at most three positional arguments as
    follows:

    - a single positional argument, if none of ``domain``,
      ``codomain``, ``values`` or ``distribution`` are specified, is
      interpreted as a FindStat map identifier.  If further arguments
      are given and it is a string, it is interpreted as a domain.
      If all this fails, it is interpreted as the specification of
      values.

    - if two positional arguments are given, the first is interpreted
      as domain, and the second either as codomain or the
      specification of values.

    - if three positional arguments are given, the first two are
      interpreted as domain and codomain, and the third as the
      specification of values.

    OUTPUT:

    An instance of a :class:`FindStatMap`, :class:`FindStatMapQuery`
    or :class:`FindStatMaps`.

    EXAMPLES:

    A particular map can be retrieved by its Mp-identifier or
    number::

        sage: findmap('Mp00062')                                                # optional -- internet
        Mp00062: Lehmer-code to major-code bijection

        sage: findmap(62)                                                       # optional -- internet
        Mp00062: Lehmer-code to major-code bijection

        sage: findmap("Mp00099oMp00127")                                        # optional -- internet
        Mp00099oMp00127

    The database can be searched by providing a list of pairs::

        sage: l = [pi for n in range(5) for pi in Permutations(n)]
        sage: findmap([(pi, pi.complement().increasing_tree_shape()) for pi in l], depth=2) # optional -- internet
        0: Mp00061oMp00069 (quality [...])

    or a dictionary::

        sage: findmap({pi: pi.complement().increasing_tree_shape() for pi in l}, depth=2) # optional -- internet
        0: Mp00061oMp00069 (quality [...])

    Note however, that the results of these two queries need not
    compare equal, because we compare queries by the data
    sent, and the ordering of the data might be different.

    Another possibility is to send a collection and a function.  In
    this case, the function is applied to the first few objects of
    the collection::

        sage: findmap("Permutations", lambda pi: pi.increasing_tree_shape(), depth=1)     # optional -- internet
        0: Mp00061 (quality [100])

    In rare cases, it may not be possible to guess the codomain of a
    map, in which case it can be provided as second argument or
    keyword argument::

        sage: findmap("Dyck paths", "Perfect matchings", lambda D: [(a+1, b) for a,b in D.tunnels()]) # optional -- internet
        0: Mp00146 (quality [100])

        sage: findmap("Dyck paths", "Set partitions", lambda D: [(a+1, b) for a,b in D.tunnels()]) # optional -- internet
        0: Mp00092oMp00146 (quality [...])

    Finally, we can also retrieve all maps with a given domain or codomain::

        sage: findmap("Cc0024")                                                 # optional -- internet
        Set of combinatorial maps with domain Cc0024: Binary words used by FindStat

        sage: findmap(codomain='Cores')                                         # optional -- internet
        Set of combinatorial maps with codomain Cc0013: Cores used by FindStat
    """
    if len(args) > 3:
        raise TypeError("findmap takes at most 3 positional arguments (%s given)" % len(args))

    bad_args = set(kwargs).difference(["values", "distribution",
                                       "domain", "codomain",
                                       "depth", "max_values"])
    if bad_args:
        raise TypeError("findmap got unexpected keyword arguments '%s'" % bad_args)

    max_values = kwargs.get("max_values", FINDSTAT_MAX_VALUES)
    depth = kwargs.get("depth", FINDSTAT_DEFAULT_DEPTH)
    values = kwargs.get("values", None)
    distribution = kwargs.get("distribution", None)
    domain = kwargs.get("domain", None)
    codomain = kwargs.get("codomain", None)

    try:
        max_values = int(max_values)
        assert 0 <= max_values <= FINDSTAT_MAX_VALUES
    except (ValueError, AssertionError):
        raise ValueError("the maximal number of values for a FindStat query must be a nonnegative integer less than or equal to %i" % FINDSTAT_MAX_VALUES)

    check_collection = True

    def get_values(raw, domain=None, codomain=None):
        if callable(raw):
            known_terms = _data_from_function(raw, domain)
            if codomain is None:
                codomain = FindStatCollection(known_terms[0][1][0])
            function = raw
        else:
            known_terms, domain, codomain = _data_from_iterable(raw, domain=domain,
                                                                codomain=codomain,
                                                                mapping=True,
                                                                check=check_collection)
            function = None
        data = _data_from_data(known_terms, max_values)
        return known_terms, data, domain, codomain, function

    def get_distribution(raw, domain=None, codomain=None):
        if callable(raw):
            known_terms = _data_from_function(raw, domain)
            function = raw
        else:
            known_terms, domain, codomain = _data_from_iterable(raw, domain=domain,
                                                             codomain=codomain,
                                                             mapping=True,
                                                             check=check_collection)
            function = None
        data = _distribution_from_data(known_terms, domain, max_values)
        return known_terms, data, domain, codomain, function

    def is_collection(arg):
        try:
            FindStatCollection(arg)
            return True
        except ValueError:
            return False

    def check_domain(arg, domain):
        if domain is not None:
            raise TypeError("the domain was specified twice, as positional argument (%s) and as keyword domain=%s" % (arg, domain))
        return arg

    def check_codomain(arg, codomain):
        if codomain is not None:
            raise TypeError("the codomain was specified twice, as positional argument (%s) and as keyword codomain=%s" % (arg, codomain))
        return arg

    def check_values(arg, values):
        if values is not None:
            raise TypeError("values were specified twice, as positional argument (%s) and as keyword values=%s" % (arg, values))
        return arg

    ######################################################################
    if values is not None and distribution is not None:
        raise ValueError("not both of `values` and `distribution` may be given for a FindStat query")

    if len(args) == 1:
        if (values is None and distribution is None
            and domain is None and codomain is None
            and (isinstance(args[0], (int, Integer, FindStatCombinatorialMap))
                 or (isinstance(args[0], str)
                     and not is_collection(args[0])))):
            return FindStatMap(args[0])

        if (isinstance(args[0], str) and
            is_collection(args[0])):
            domain = check_domain(args[0], domain)

        else:
            values = check_values(args[0], values)

    elif len(args) == 2:
        domain = check_domain(args[0], domain)
        if isinstance(args[1], (int, Integer, str)):
            codomain = check_codomain(args[1], codomain)
        else:
            values = check_values(args[1], values)

    elif len(args) == 3:
        domain = check_domain(args[0], domain)
        codomain = check_codomain(args[1], codomain)
        values = check_values(args[2], values)

    if domain is not None:
        domain = FindStatCollection(domain)
    if codomain is not None:
        codomain = FindStatCollection(codomain)

    if (values is None and distribution is None
        and (domain is not None or codomain is not None)):
        return FindStatMaps(domain=domain, codomain=codomain)

    if values is not None:
        if isinstance(values, (int, Integer, str, FindStatCombinatorialMap)):
            if domain is not None or codomain is not None:
                raise ValueError("domain and codomain must not be provided if a map identifier is given")
            return FindStatMapQuery(values_of=values, depth=depth)

        known_terms, data, domain, codomain, function = get_values(values, domain, codomain)
        return FindStatMapQuery(data=data, domain=domain, codomain=codomain, depth=depth,
                                known_terms=known_terms, function=function)

    if distribution is not None:
        if isinstance(distribution, (int, Integer, str, FindStatCombinatorialMap)):
            if domain is not None or codomain is not None:
                raise ValueError("domain and codomain must not be provided if a map identifier is given")
            return FindStatMapQuery(distribution_of=distribution, depth=depth)

        known_terms, data, domain, function = get_distribution(distribution, domain)
        return FindStatMapQuery(data=data, domain=domain, codomain=codomain, depth=depth,
                                known_terms=known_terms, function=function)

    raise ValueError("the given arguments cannot be used for a FindStat search")


######################################################################
# common methods for statistic and maps in FindStat
######################################################################
class FindStatFunction(SageObject):
    """
    A class providing the common methods of :class:`FindStatMap` and
    :class:`FindStatStatistic`.

    This class provides methods to access and modify properties of a
    single statistic or map of the FindStat database.
    """
    def __init__(self, id, data=None, function=None):
        """
        Initialize a statistic or map.

        INPUT:

        - ``id`` -- a padded identifier, with number 0 reserved for new
          statistics or maps.

        - ``data`` -- dictionary with "Description", "Code", etc.

        - ``function`` -- (optional) a callable implementing the
          statistic or map, or ``None``

        ``data`` should be provided if and only if ``id`` refers to a
        new statistic or map (with identifier 0).

        TESTS::

            sage: from sage.databases.findstat import FindStatFunction, FindStatCollection
            sage: FindStatFunction("St000000",                                  # optional -- internet
            ....:                  data={"Bibliography": {},
            ....:                        "Code": "",
            ....:                        "Description" : "",
            ....:                        "Domain": FindStatCollection(1),
            ....:                        "Name": "a new statistic",
            ....:                        "References": "",
            ....:                        "SageCode": ""})
            St000000: a new statistic
        """
        self._id = id # as padded identifier, with number 0 reserved for new statistics or maps
        self._modified = False # set in every method modifying the data
        if callable(function):
            self._function = function
        else:
            self._function = False  # determines that FindStat code may not be executed
        if self.id() != 0 and data is not None:
            raise ValueError("data (%s) may be provided if and only if id (%s) is %s or %s" %
                             (data, id,
                              FINDSTAT_STATISTIC_PADDED_IDENTIFIER % 0,
                              FINDSTAT_MAP_PADDED_IDENTIFIER % 0))
        self._data_cache = data # a dictionary with "Description", "Code", etc.

    def _data(self):
        """
        Return a copy of the data defining the statistic or map.

        The first terms of a statistic are provided separately, by
        :meth:`first_terms_raw`, to save bandwidth.

        Any method modifying ``self._data_cache`` should first
        call - directly or indirectly - this method.

        TESTS::

            sage: from sage.databases.findstat import FindStatFunction, FindStatCollection
            sage: FindStatFunction("St000000",                                  # optional -- internet
            ....:                  data={"Bibliography": {},
            ....:                        "Code": "",
            ....:                        "Description" : "",
            ....:                        "Domain": FindStatCollection(1),
            ....:                        "Name": "a new statistic",
            ....:                        "References": "",
            ....:                        "SageCode": ""})._data()
            {'Bibliography': {},
             'Code': '',
             'Description': '',
             'Domain': Cc0001: Permutations,
             'Name': 'a new statistic',
             'References': '',
             'SageCode': ''}
        """
        # initializes self._data_cache on first call
        if self._data_cache is None:
            self._data_cache = self._fetch_data()
        # some of the data are lists, so we need to deepcopy
        return deepcopy(self._data_cache)

    def __call__(self, elt):
        """
        Apply the function to a given element.

        EXAMPLES::

            sage: s = lambda g: g.diameter() if g.is_connected() else None
            sage: q = findstat("graphs", s, max_values=100)                     # optional -- internet
            sage: q(graphs.PetersenGraph().copy(immutable=True))                # optional -- internet
            2
        """
        if self._function is False and FindStat()._allow_execution is False:
            raise ValueError("execution of verified code provided by FindStat is not enabled for %s" % self)
        if self._function is True or (self._function is False and FindStat()._allow_execution is True):
            if not self.sage_code():
                raise ValueError("there is no verified code available for %s" % self)
            from sage.repl.preparse import preparse
            try:
                l = {}
                environment = 'sage.all'
                code = f"from {environment} import *\n" + preparse(self.sage_code())
                exec(code, l)
            except SyntaxError:
                raise ValueError("could not execute verified code for %s" % self)
            if isinstance(self, FindStatStatistic):
                self._function = eval("statistic", l)
            elif isinstance(self, FindStatMap):
                self._function = eval("mapping", l)
            else:
                raise ValueError("cannot execute verified code for %s" % self)

        return self._function(elt)

    def __repr__(self):
        r"""
        Return the representation of the FindStat statistic or map.

        OUTPUT:

        A string, the identifier and the name of the statistic.  If
        the statistic was modified (see :meth:`modified`) this is
        also indicated.

        EXAMPLES::

            sage: findstat(51)                                                  # optional -- internet
            St000051: The size of the left subtree of a binary tree.

            sage: findstat(914)                                                 # optional -- internet
            St000914: The sum of the values of the Möbius function of a poset.

            sage: findmap(85)                                                   # optional -- internet
            Mp00085: Schützenberger involution
        """
        if self._modified:
            s = "%s(modified): %s" % (self.id_str(), self.name())
        else:
            s = "%s: %s" % (self.id_str(), self.name())
        return s

    def reset(self):
        """
        Discard all modification of the statistic or map.

        EXAMPLES::

            sage: s = findmap(62)                                               # optional -- internet
            sage: s.set_name(u"Möbius"); s                                      # optional -- internet
            Mp00062(modified): Möbius
            sage: s.reset(); s                                                  # optional -- internet
            Mp00062: Lehmer-code to major-code bijection

        TESTS:

        Check that new statistics and maps cannot be reset::

            sage: # optional - internet
            sage: q = findstat([(d, randint(1, 1000)) for d in DyckWords(4)])
            sage: q.set_description("Random values on Dyck paths.")
            sage: print(q.description())
            Random values on Dyck paths.
            sage: q.reset()
            Traceback (most recent call last):
            ...
            ValueError: cannot reset values of St000000: a new statistic on Dyck paths
        """
        if isinstance(self, FindStatStatistic) and self.id_str() in _all_statistics:
            del _all_statistics[self.id_str()]
        elif isinstance(self, FindStatMap) and self.id_str() in _all_maps:
            del _all_maps[self.id_str()]
        else:
            raise ValueError("cannot reset values of %s" % self)
        self.__init__(self.parent(), self.id_str())

    def id(self):
        r"""
        Return the FindStat identifier of the statistic or map.

        OUTPUT: the FindStat identifier of the statistic or map, as an integer

        EXAMPLES::

            sage: findstat(51).id()                                             # optional -- internet
            51
        """
        return int(self._id[2:])

    def id_str(self):
        r"""
        Return the FindStat identifier of the statistic or map.

        OUTPUT: the FindStat identifier of the statistic or map, as a string

        EXAMPLES::

            sage: findstat(51).id_str()                                         # optional -- internet
            'St000051'
        """
        return self._id

    def domain(self):
        r"""
        Return the FindStat domain of the statistic or map.

        OUTPUT:

        The domain of the statistic or map as an instance of
        :class:`FindStatCollection`.

        EXAMPLES::

            sage: findstat(51).domain()                                         # optional -- internet
            Cc0010: Binary trees

            sage: findmap(62).domain()                                          # optional -- internet
            Cc0001: Permutations
        """
        return FindStatCollection(self._data()["Domain"])

    def description(self):
        r"""
        Return the description of the statistic or map.

        OUTPUT: string; for statistics, the first line is used as name

        EXAMPLES::

            sage: print(findstat(51).description())                             # optional -- internet
            The size of the left subtree of a binary tree.
        """
        return self._data()["Description"]

    def set_description(self, value):
        r"""
        Set the description of the statistic or map.

        INPUT:

        - ``value`` -- string; for statistics, this is the name of the
          statistic followed by its description on a separate line

        This information is used when submitting the statistic or map with
        :meth:`submit`.

        EXAMPLES::

            sage: q = findstat([(d, randint(1, 1000)) for d in DyckWords(4)])             # optional -- internet
            sage: q.set_description("Random values on Dyck paths.\nNot for submission.")  # optional -- internet
            sage: print(q.description())                                                  # optional -- internet
            Random values on Dyck paths.
            Not for submission.
        """
        if value != self.description():
            self._modified = True
            self._data_cache["Description"] = value

    def name(self):
        r"""
        Return the name of the statistic or map.

        OUTPUT:

        A string.  For statistics, this is just the first line of the
        description.

        EXAMPLES::

            sage: findstat(51).name()                                           # optional -- internet
            'The size of the left subtree of a binary tree.'
        """
        return self._data()["Name"]

    def references(self):
        r"""
        Return the references associated with the statistic or map.

        OUTPUT:

        An instance of :class:`sage.databases.oeis.FancyTuple`, each
        item corresponds to a reference.

        EXAMPLES::

            sage: findstat(41).references()                                     # optional -- internet
            0: [1]  de Médicis, A., Viennot, X. G., Moments des $q$-polynômes de Laguerre et la bijection de Foata-Zeilberger [[MathSciNet:1288802]]
            1: [2]  Simion, R., Stanton, D., Octabasic Laguerre polynomials and permutation statistics [[MathSciNet:1418763]]
        """
        result = []
        refs = self.references_raw()
        if refs:
            refs = refs.splitlines()
        else:
            return FancyTuple([])
        bibs = self._data()["Bibliography"]
        for ref in refs:
            parts = ref.partition("[[")
            parts = parts[:-1] + parts[2].partition("]]")
            comment = parts[0]
            link = parts[2]
            if link == "":
                result.append(ref)
            else:
                try:
                    bibitem = bibs[link]
                except KeyError:
                    # this means that the link is unhandled
                    result.append(ref)
                else:
                    author_title = ", ".join(e for e in [bibitem["Author"], bibitem["Title"]]
                                             if e)
                    result.append(comment + author_title + " " + "".join(parts[1:]))

        return FancyTuple(result)

    def references_raw(self):
        r"""
        Return the unrendered references associated with the statistic or map.

        EXAMPLES::

            sage: print(findstat(41).references_raw())                          # optional -- internet
            [1]  [[MathSciNet:1288802]]
            [2]  [[MathSciNet:1418763]]
        """
        return self._data()["References"]

    def set_references_raw(self, value):
        r"""
        Set the references associated with the statistic or map.

        INPUT:

        - ``value`` -- string; each reference should be on a single line, and
          consist of one or more links to the same item

        FindStat will automatically resolve the links, if possible.
        A complete list of supported services can be found at
        <https://findstat.org/NewStatistic>.

        This information is used when submitting the statistic with
        :meth:`submit`.

        EXAMPLES::

            sage: q = findstat([(d, randint(1, 1000)) for d in DyckWords(4)])   # optional -- internet
            sage: q.set_references_raw("[[arXiv:1102.4226]]\n[[oeis:A000001]]") # optional -- internet
            sage: q.references()                                                # optional -- internet
            0: [[arXiv:1102.4226]]
            1: [[oeis:A000001]]
        """
        if value != self.references_raw():
            self._modified = True
            self._data_cache["References"] = value

    def sage_code(self):
        r"""
        Return the Sage code associated with the statistic or map.

        OUTPUT: an empty string or a string of the form::

            def statistic(x):
                ...

        or::

            def mapping(x):
                ...

        EXAMPLES::

            sage: print(findstat(51).sage_code())                               # optional -- internet
            def statistic(T):
                return T[0].node_number()
        """
        return self._data()["SageCode"]

    def set_sage_code(self, value):
        r"""
        Set the code associated with the statistic or map.

        INPUT:

        - ``value`` -- string; SageMath code producing the values of the
          statistic or map

        Contributors are encouraged to submit code for statistics
        using :meth:`FindStatStatistic.set_code`.  Modifying the
        "verified" SageMath code using this method is restricted to
        members of the FindStatCrew, for all other contributors this
        method has no effect.

        EXAMPLES::

            sage: q = findstat([(d, randint(1,1000)) for d in DyckWords(4)])              # optional -- internet
            sage: q.set_sage_code("def statistic(x):\n    return randint(1,1000)")        # optional -- internet
            sage: print(q.sage_code())                                                    # optional -- internet
            def statistic(x):
                return randint(1,1000)
        """
        if value != self.sage_code():
            self._modified = True
            self._data_cache["SageCode"] = value

######################################################################
# statistics
######################################################################


class FindStatCombinatorialStatistic(SageObject):
    """
    A class providing methods to retrieve the first terms of a statistic.

    This class provides methods applicable to instances of
    :class:`FindStatStatistic`, :class:`FindStatCompoundStatistic`
    and :class:`FindStatStatisticQuery`.
    """
    def __init__(self):
        """
        Initialize the combinatorial statistic.

        TESTS::

            sage: from sage.databases.findstat import FindStatCombinatorialStatistic
            sage: FindStatCombinatorialStatistic()
            <sage.databases.findstat.FindStatCombinatorialStatistic object at 0x...>
        """
        self._first_terms_cache = None
        self._first_terms_raw_cache = None

    def first_terms(self):
        r"""
        Return the first terms of the (compound) statistic as a
        dictionary.

        OUTPUT:

        A dictionary from Sage objects representing an element of the
        appropriate collection to integers.

        This method is overridden in :class:`FindStatStatisticQuery`.

        EXAMPLES::

            sage: findstat(41).first_terms()[PerfectMatching([(1,6),(2,5),(3,4)])]        # optional -- internet
            3
        """
        # initialize self._first_terms_cache and
        # self._first_terms_raw_cache on first call
        if self._first_terms_cache is None:
            self._first_terms_cache = self._fetch_first_terms()
        # a shallow copy suffices - tuples are immutable
        return dict(self._first_terms_cache)

    def _first_terms_raw(self, max_values):
        """
        Return the first terms of the (compound) statistic as a
        list of ``(string, value)`` pairs.

        INPUT:

        - ``max_values`` -- integer determining how many terms to
          return at most

        OUTPUT: list of ``(string, value)`` pairs

        This method is overridden in :class:`FindStatStatisticQuery`.

        TESTS::

            sage: findstat(41)._first_terms_raw(4)                              # optional -- internet
            [('[(1,2)]', 0),
             ('[(1,2),(3,4)]', 0),
             ('[(1,3),(2,4)]', 0),
             ('[(1,4),(2,3)]', 1)]
        """
        # initialize self._first_terms_raw_cache on first call
        if self._first_terms_raw_cache is None:
            self._first_terms_raw_cache = self._fetch_first_terms_raw()
        # a shallow copy suffices - tuples are immutable
        return self._first_terms_raw_cache[:max_values]

    def first_terms_str(self, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Return the first terms of the statistic in the format needed
        for a FindStat query.

        OUTPUT:

        A string, where each line is of the form ``object => value``,
        where ``object`` is the string representation of an element
        of the appropriate collection as used by FindStat and value
        is an integer.

        EXAMPLES::

            sage: print(findstat(41).first_terms_str(max_values=4))             # optional -- internet
            [(1,2)] => 0
            [(1,2),(3,4)] => 0
            [(1,3),(2,4)] => 0
            [(1,4),(2,3)] => 1

        TESTS:

        Check that no more terms than asked for are computed::

            sage: st = cached_function(lambda d: randint(1,1000))
            sage: s = findstat("Dyck paths", st, max_values=100, depth=0)       # optional -- internet
            sage: len(st.cache)                                                 # optional -- internet
            100
            sage: _ = s.first_terms_str(max_values=100)                         # optional -- internet
            sage: len(st.cache)                                                 # optional -- internet
            100
        """
        return "\n".join(key + " => " + str(val)
                         for key, val in self._first_terms_raw(max_values=max_values))

    def _fetch_first_terms(self):
        r"""
        Return the first terms of the statistic as a list of ``(object,
        value)`` pairs.

        TESTS::

            sage: findstat(41)._fetch_first_terms()[:4]                         # optional -- internet
            [([(1, 2)], 0),
             ([(1, 2), (3, 4)], 0),
             ([(1, 3), (2, 4)], 0),
             ([(1, 4), (2, 3)], 1)]
        """
        from_str = self.domain().from_string()
        if self._first_terms_raw_cache is None:
            self._first_terms_raw_cache = self._fetch_first_terms_raw()
        return [(from_str(obj), Integer(val))
                for obj, val in self._first_terms_raw_cache]

    def _generating_functions_dict(self,
                                   max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Return the generating functions of ``self`` as dictionary of
        dictionaries, computed from ``self.first_terms``.

        TESTS:

            sage: q = findstat((BinaryTrees(5), list(range(0,24,4))*7))         # optional -- internet
            sage: q._generating_functions_dict()                                # optional -- internet
            {5: {0: 7, 4: 7, 8: 7, 12: 7, 16: 7, 20: 7}}
        """
        gfs = {}
        lvls = {}
        domain = self.domain()
        levels_with_sizes = domain.levels_with_sizes()
        total = 0
        for elt, val in self.first_terms().items():
            if total == max_values:
                break
            lvl = domain.element_level(elt)
            if lvl not in levels_with_sizes:
                continue
            total += 1
            if lvl not in gfs:
                gfs[lvl] = {}
            gfs[lvl][val] = gfs[lvl].get(val, 0) + 1
            lvls[lvl] = lvls.get(lvl, 0) + 1

        for lvl, size in lvls.items():
            if size < levels_with_sizes[lvl]:
                del gfs[lvl]
        return gfs

    def generating_functions(self, style='polynomial',
                             max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Return the generating functions of the statistic as a dictionary.

        The keys of this dictionary are the levels for which the
        generating function of the statistic can be computed from
        the known data.  Each value represents a generating function
        for one level, as a polynomial, as a dictionary, or as a list
        of coefficients.

        INPUT:

        - ``style`` -- string (default: ``'polynomial'``); can be
          ``'polynomial'``, ``'dictionary'``, or ``'list'``

        OUTPUT:

        - if ``style`` is ``'polynomial'``, the generating function is
          returned as a polynomial

        - if ``style`` is ``'dictionary'``, the generating function is
          returned as a dictionary representing the monomials of the
          generating function

        - if ``style`` is ``'list'``, the generating function is
          returned as a list of coefficients of the generating
          function.  In this case, leading and trailing zeros are
          omitted.

        EXAMPLES::

            sage: st = findstat(41)                                             # optional -- internet
            sage: st.generating_functions()                                     # optional -- internet
            {2: 1,
             4: q + 2,
             6: q^3 + 3*q^2 + 6*q + 5,
             8: q^6 + 4*q^5 + 10*q^4 + 20*q^3 + 28*q^2 + 28*q + 14}

            sage: st.generating_functions(style='dictionary')                   # optional -- internet
            {2: {0: 1},
             4: {0: 2, 1: 1},
             6: {0: 5, 1: 6, 2: 3, 3: 1},
             8: {0: 14, 1: 28, 2: 28, 3: 20, 4: 10, 5: 4, 6: 1}}

            sage: st.generating_functions(style='list')                         # optional -- internet
            {2: [1], 4: [2, 1], 6: [5, 6, 3, 1], 8: [14, 28, 28, 20, 10, 4, 1]}

        TESTS::

            sage: st = findstat(41)                                             # optional -- internet
            sage: st.generating_functions(max_values=19)                        # optional -- internet
            {2: 1, 4: q + 2, 6: q^3 + 3*q^2 + 6*q + 5}

            sage: st = findstat("graphs", lambda G: G.size(), max_values=100)   # optional -- internet
            sage: st.generating_functions(max_values=18)                        # optional -- internet
            {1: 1,
             2: q + 1,
             3: q^3 + q^2 + q + 1,
             4: q^6 + q^5 + 2*q^4 + 3*q^3 + 2*q^2 + q + 1}
            sage: st.generating_functions(max_values=1252)                      # optional -- internet
            {1: 1,
             2: q + 1,
             3: q^3 + q^2 + q + 1,
             4: q^6 + q^5 + 2*q^4 + 3*q^3 + 2*q^2 + q + 1,
             5: q^10 + q^9 + 2*q^8 + 4*q^7 + 6*q^6 + 6*q^5 + 6*q^4 + 4*q^3 + 2*q^2 + q + 1,
             6: q^15 + q^14 + 2*q^13 + 5*q^12 + 9*q^11 + 15*q^10 + 21*q^9 + 24*q^8 + 24*q^7 + 21*q^6 + 15*q^5 + 9*q^4 + 5*q^3 + 2*q^2 + q + 1,
             7: q^21 + q^20 + 2*q^19 + 5*q^18 + 10*q^17 + 21*q^16 + 41*q^15 + 65*q^14 + 97*q^13 + 131*q^12 + 148*q^11 + 148*q^10 + 131*q^9 + 97*q^8 + 65*q^7 + 41*q^6 + 21*q^5 + 10*q^4 + 5*q^3 + 2*q^2 + q + 1}
        """
        d = self._generating_functions_dict(max_values=max_values)
        return _generating_functions_from_dict(d, style)

    def oeis_search(self, search_size=32, verbose=True):
        r"""
        Search the OEIS for the generating function of the statistic.

        INPUT:

        - ``search_size`` -- (default: 32) the number of integers in the
          sequence. If this is chosen too big, the OEIS result may be
          corrupted.

        - ``verbose`` -- boolean (default: ``True``); if ``True``, some
          information about the search are printed

        OUTPUT: a tuple of OEIS sequences, see
        :meth:`sage.databases.oeis.OEIS.find_by_description` for more
        information

        EXAMPLES::

            sage: st = findstat(41)                                             # optional -- internet

            sage: st.oeis_search()                                              # optional -- internet
            Searching the OEIS for "1  2,1  5,6,3,1  14,28,28,20,10,4,1"
            0: A067311: Triangle read by rows: T(n,k) gives number of ways of arranging n chords on a circle with k simple intersections ...
        """
        from sage.databases.oeis import oeis
        gen_funcs = self.generating_functions(style='list')

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
            OEIS_func_string = ",".join(str(coefficient) for coefficient in gen_func)
            OEIS_string += OEIS_func_string + "  "
        OEIS_string = OEIS_string.strip()
        if counter >= 4:
            if verbose:
                print('Searching the OEIS for "%s"' % OEIS_string)
            return oeis(OEIS_string)

        if verbose:
            print("Too little information to search the OEIS for this statistic (only %s values given)." % counter)


class FindStatStatistic(Element,
                        FindStatFunction,
                        FindStatCombinatorialStatistic,
                        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A FindStat statistic.

    :class:`FindStatStatistic` is a class representing a
    combinatorial statistic available in the FindStat database.

    This class provides methods to inspect and update various
    properties of these statistics.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatStatistic
        sage: FindStatStatistic(41)                                             # optional -- internet
        St000041: The number of nestings of a perfect matching.

    .. SEEALSO::

        :class:`FindStatStatistics`
    """
    @staticmethod
    def __classcall_private__(cls, entry):
        """
        Retrieve a statistic from the database.

        TESTS::

            sage: from sage.databases.findstat import FindStatStatistic
            sage: FindStatStatistic("abcdefgh")                                 # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: the value 'abcdefgh' is not a valid FindStat statistic identifier
        """
        return FindStatStatistics()(entry)

    def __init__(self, parent, id):
        """
        Initialize a FindStat statistic from an identifier.

        INPUT:

        - ``parent`` -- :class:`FindStatStatistics`

        - ``id`` -- the (padded) FindStat identifier of the statistic

        EXAMPLES::

            sage: findstat(41)                                                  # optional -- internet, indirect doctest
            St000041: The number of nestings of a perfect matching.
        """
        FindStatFunction.__init__(self, id)
        FindStatCombinatorialStatistic.__init__(self)
        Element.__init__(self, parent)

    def __call__(self, elt):
        """
        Apply the statistic to a given element.

        EXAMPLES::

            sage: s = lambda g: g.diameter() if g.is_connected() else None
            sage: q = findstat("graphs", s, max_values=100)                     # optional -- internet
            sage: q(graphs.PetersenGraph().copy(immutable=True))                # optional -- internet
            2
        """
        val = self.first_terms().get(elt, None)
        if val is None:
            return FindStatFunction.__call__(self, elt)
        return val

    def __reduce__(self):
        """
        Return a function and its arguments needed to create the
        statistic.

        TESTS::

            sage: from sage.databases.findstat import FindStatStatistic
            sage: c = FindStatStatistic(41)                                     # optional -- internet
            sage: loads(dumps(c)) == c                                          # optional -- internet
            True
        """
        return FindStatStatistic, (self.id(),)

    def _richcmp_(self, other, op):
        """
        Compare two statistics by identifier.

        TESTS::

            sage: findstat(41) != findstat(42)                                  # optional -- internet
            True
            sage: findstat(41) == findstat(41)                                  # optional -- internet
            True
        """
        return richcmp(self.id(), other.id(), op)

    def _fetch_data(self):
        r"""
        Return a dictionary containing the data of the statistic, except
        for the values, fetched from FindStat.

        TESTS::

            sage: findstat(41)._data()                                          # optional -- internet, indirect doctest
            {'Bibliography': {'MathSciNet:1288802': {'Author': 'de Médicis, A., Viennot, X. G.',
               'Title': 'Moments des $q$-polynômes de Laguerre et la bijection de Foata-Zeilberger'},
              'MathSciNet:1418763': {'Author': 'Simion, R., Stanton, D.',
               'Title': 'Octabasic Laguerre polynomials and permutation statistics'}},
             'Code': 'def statistic(x):\r\n    return len(x.nestings())',
             'Description': 'The number of nestings of a perfect matching.\r\n\r\nThis is the number of pairs of edges $((a,b), (c,d))$ such that $a\\le c\\le d\\le b$. i.e., the edge $(c,d)$ is nested inside $(a,b)$.',
             'Domain': 'Cc0012',
             'Name': 'The number of nestings of a perfect matching.',
             'References': '[1]  [[MathSciNet:1288802]]\r\n[2]  [[MathSciNet:1418763]]',
             'SageCode': 'def statistic(x):\r\n    return len(x.nestings())'}
        """
        fields = "Bibliography,Code,Description,Domain,Name,References,SageCode"
        fields_Bibliography = "Author,Title"
        url = (FINDSTAT_API_STATISTICS + self.id_str()
               + "?fields=" + fields
               + "&fields[Bibliography]=" + fields_Bibliography)
        verbose("fetching statistic data %s" % url, caller_name='FindStatStatistic')

        included = _get_json(url)["included"]
        # slightly simplify the representation
        data = {key: val for key, val in included["Statistics"][self.id_str()].items()}
        # we replace the list of identifiers in Bibliography with the dictionary
        data["Bibliography"] = included["References"]
        return data

    def _fetch_first_terms_raw(self):
        r"""
        Return the first terms of the statistic, as ``(string,
        value)`` pairs, fetched from FindStat

        TESTS::

            sage: findstat(41)._first_terms_raw(4)                              # optional -- internet, indirect doctest
            [('[(1,2)]', 0),
             ('[(1,2),(3,4)]', 0),
             ('[(1,3),(2,4)]', 0),
             ('[(1,4),(2,3)]', 1)]
        """
        fields = "Values"
        url = FINDSTAT_API_STATISTICS + self.id_str() + "?fields=" + fields
        values = _get_json(url)["included"]["Statistics"][self.id_str()]["Values"]
        return [tuple(pair) for pair in values]

    def set_first_terms(self, values):
        r"""
        Update the first terms of the statistic.

        INPUT:

        - ``values`` -- list of pairs of the form ``(object, value)`` where
          ``object`` is a Sage object representing an element of the
          appropriate collection and ``value`` is an integer

        This information is used when submitting the statistic with
        :meth:`submit`.

        .. WARNING::

            This method cannot check whether the given values are
            actually correct.  Moreover, it does not even perform any
            sanity checks.

        TESTS::

            sage: s = findstat(41)                                              # optional -- internet
            sage: l = [([(1,2)], 1), ([(1,2),(3,4)], 7), ([(1,3),(2,4)], 8), ([(1,4),(2,3)], 3)]
            sage: s.set_first_terms(l)                                          # optional -- internet
            sage: print(s.first_terms_str())                                    # optional -- internet
            [(1, 2)] => 1
            [(1, 2), (3, 4)] => 7
            [(1, 3), (2, 4)] => 8
            [(1, 4), (2, 3)] => 3
            sage: s.reset()                                                     # optional -- internet
        """
        to_str = self.domain().to_string()
        new = [(to_str(obj), value) for obj, value in values]
        if sorted(new) != sorted(self.first_terms_str()):
            self._modified = True
            self._first_terms_raw_cache = new
            self._first_terms_cache = values

    def code(self):
        r"""
        Return the code associated with the statistic or map.

        OUTPUT: string

        Contributors are encouraged to submit Sage code in the form::

            def statistic(x):
                ...

        but the string may also contain code for other computer
        algebra systems.

        EXAMPLES::

            sage: print(findstat(41).code())                                    # optional -- internet
            def statistic(x):
                return len(x.nestings())

            sage: print(findstat(118).code())                                   # optional -- internet
            (* in Mathematica *)
            tree = {{{{}, {}}, {{}, {}}}, {{{}, {}}, {{}, {}}}};
            Count[tree, {{___}, {{___}, {{___}, {___}}}}, {0, Infinity}]
        """
        return self._data()["Code"]

    def set_code(self, value):
        r"""
        Set the code associated with the statistic.

        INPUT:

        - ``value`` -- string; code producing the values of the statistic

        Contributors are encouraged to submit SageMath code in the form::

            def statistic(x):
                ...

        However, code for any other platform is accepted also.

        This information is used when submitting the statistic with
        :meth:`submit`.

        EXAMPLES::

            sage: q = findstat([(d, randint(1,1000)) for d in DyckWords(4)])    # optional -- internet
            sage: q.set_code("def statistic(x):\n    return randint(1,1000)")   # optional -- internet
            sage: print(q.code())                                               # optional -- internet
            def statistic(x):
                return randint(1,1000)
        """
        if value != self.code():
            self._modified = True
            self._data_cache["Code"] = value

    def browse(self):
        r"""
        Open the FindStat web page of the statistic in a browser.

        EXAMPLES::

            sage: findstat(41).browse()                                         # optional -- webbrowser
        """
        if self.id() == 0:
            self.submit()
        else:
            webbrowser.open(FINDSTAT_URL_STATISTICS + self.id_str())

    def submit(self, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        r"""
        Open the FindStat web page for editing the statistic or
        submitting a new statistic in a browser.

        TESTS::

            sage: s = findstat([(d, randint(1,1000)) for d in DyckWords(4)])    # optional -- internet
            sage: s.set_description(u"Möbius")                                  # optional -- internet
            sage: s.submit()                                                    # optional -- webbrowser
        """
        args = {}
        args["OriginalStatistic"] = self.id_str()
        args["Domain"] = self.domain().id_str()
        args["Values"] = self.first_terms_str(max_values=max_values)
        args["Description"] = self.description()
        args["References"] = self.references_raw()
        args["Code"] = self.code()
        args["SageCode"] = self.sage_code()
        args["CurrentAuthor"] = FindStat().user_name()
        args["CurrentEmail"] = FindStat().user_email()

        if not self.id():
            url = FINDSTAT_NEWSTATISTIC_FORM_HEADER % FINDSTAT_URL_NEW_STATISTIC
        else:
            url = FINDSTAT_NEWSTATISTIC_FORM_HEADER % (FINDSTAT_URL_EDIT_STATISTIC + self.id_str())
        _submit(args, url)

    # editing and submitting is really the same thing
    edit = submit

    def __hash__(self):
        """
        Return a hash value for the statistic.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: list(FindStatMaps(domain=1, codomain=10))                     # optional -- internet
            [Mp00061: to increasing tree, Mp00072: binary search tree: left to right]
        """
        return self.id()

    def info(self):
        """
        Print a detailed description of the statistic.

        EXAMPLES::

            sage: findstat("St000042").info()                                   # optional -- internet
                St000042: The number of crossings of a perfect matching.
        """
        print("    %s" % self)


_all_statistics = {}


class FindStatStatistics(UniqueRepresentation, Parent):
    r"""
    The class of FindStat statistics.

    The elements of this class are combinatorial statistics currently
    in FindStat.

    EXAMPLES:

    We can print a list of the first few statistics currently in
    FindStat in a given domain::

        sage: from sage.databases.findstat import FindStatStatistics
        sage: for st, _ in zip(FindStatStatistics("Perfect Matchings"), range(3)):        # optional -- internet
        ....:     print("    " + st.name())
        The number of nestings of a perfect matching.
        The number of crossings of a perfect matching.
        The number of crossings plus two-nestings of a perfect matching.
    """
    def __init__(self, domain=None):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatStatistics
            sage: M = FindStatStatistics()                                      # optional -- internet
            sage: TestSuite(M).run()                                            # optional -- internet
        """
        if domain is None:
            self._domain = None
        else:
            self._domain = FindStatCollection(domain)
        self._identifiers = None
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, id):
        """
        Initialize a FindStat statistic.

        INPUT:

        - ``id`` -- string containing the FindStat identifier of
          the statistic, or the corresponding integer

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatStatistic
            sage: FindStatStatistic(41)                                         # optional -- internet
            St000041: The number of nestings of a perfect matching.
        """
        if isinstance(id, self.Element):
            return id
        if isinstance(id, (int, Integer)):
            id = FINDSTAT_STATISTIC_PADDED_IDENTIFIER % id
        elif isinstance(id, FindStatCombinatorialStatistic):
            id = id.id_str()
        if not isinstance(id, str):
            raise TypeError("the value '%s' is not a valid FindStat statistic identifier, nor a FindStat statistic query" % id)
        else:
            id = id.strip()
        if FINDSTAT_MAP_SEPARATOR in id:
            return FindStatCompoundStatistic(id)
        if not re.match(FINDSTAT_STATISTIC_REGEXP, id) or int(id[2:]) <= 0:
            raise ValueError("the value '%s' is not a valid FindStat statistic identifier" % id)
        if id not in _all_statistics or _all_statistics[id] is None:
            _all_statistics[id] = self.element_class(self, id)

        return _all_statistics[id]

    def _repr_(self):
        """
        Return a short description of the set of FindStat statistics.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatStatistics
            sage: FindStatStatistics()                                          # optional -- internet
            Set of combinatorial statistics in FindStat

            sage: FindStatStatistics(12)                                        # optional -- internet
            Set of combinatorial statistics with domain Cc0012: Perfect matchings in FindStat
        """
        if self._domain is None:
            return "Set of combinatorial statistics in FindStat"
        return "Set of combinatorial statistics with domain %s in FindStat" % self._domain

    def __iter__(self):
        """
        Return an iterator over all FindStat statistics.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatStatistics
            sage: next(iter(FindStatStatistics("Perfect Matchings")))           # optional -- internet
            St000041: The number of nestings of a perfect matching.
        """
        if self._identifiers is None:
            if self._domain is None:
                url = FINDSTAT_API_STATISTICS
            else:
                url = FINDSTAT_API_STATISTICS + "?Domain=%s" % self._domain.id_str()

            self._identifiers = _get_json(url)["data"]

        for st in self._identifiers:
            yield FindStatStatistic(st)

    def _an_element_(self):
        """
        Return a FindStat statistic.

        EXAMPLES::

            sage: findstat(domain='Permutations').an_element()                  # optional -- internet
            St000001: The number of reduced words for a permutation.
        """
        try:
            return next(iter(self))
        except StopIteration:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError

    Element = FindStatStatistic


class FindStatStatisticQuery(FindStatStatistic):
    """
    A class representing a query for FindStat (compound) statistics.
    """
    def __init__(self, data=None, values_of=None, distribution_of=None,
                 domain=None, known_terms=None, function=None,
                 depth=FINDSTAT_DEFAULT_DEPTH,
                 debug=False):
        """
        Initialize a query for FindStat (compound) statistics.

        INPUT:

        - ``data`` -- (optional), a list of pairs ``(objects,
          values)``, where ``objects`` and ``values`` are all lists
          of the same length, the former are elements in the FindStat
          collection, the latter are integers

        - ``known_terms`` -- (optional), a lazy list in the same format
          as ``data``, which agrees with ``data``, and may be used
          for submission

        - ``values_of`` -- (optional), anything accepted by
          :class:`FindStatCompoundStatistic`

        - ``distribution_of`` -- (optional), anything accepted by
          :class:`FindStatCompoundStatistic`

        - ``domain`` -- (optional), anything accepted by
          :class:`FindStatCollection`

        - ``depth`` -- (optional), the number of maps to apply before
          applying the statistic


        Only one of ``data``, ``values_of`` and ``distribution_of``
        may be provided.  The parameter ``domain`` must be provided
        if and only if ``data`` is provided, or ``values_of`` or
        ``distribution_of`` are given as a function.

        The parameter ``known_terms`` is only allowed, if ``data`` is
        provided.  It defaults to ``data``.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatStatisticQuery
            sage: data = [[[m], [m.number_of_nestings()]] for n in range(5) for m in PerfectMatchings(2*n)]
            sage: FindStatStatisticQuery(domain=12, data=data, depth=1)         # optional -- internet
            0: St000041 (quality [99, 100])
            1: St000041oMp00113 (quality [99, 100])
            2: St000042oMp00116 (quality [99, 100])
            ...
        """
        self._first_terms = data
        if data is not None and known_terms is None:
            self._known_terms = data
        else:
            self._known_terms = known_terms
        self._values_of = None
        self._distribution_of = None
        self._depth = depth

        if data is not None:
            assert all(param is None for param in [distribution_of, values_of])

            domain = FindStatCollection(domain)
            query = {"Domain": domain.id_str(),
                     "Data": _data_to_str(self._first_terms, domain)}

        elif distribution_of is not None:
            assert all(param is None for param in [data, known_terms, values_of])

            self._distribution_of = FindStatCompoundStatistic(distribution_of)
            domain = self._distribution_of.domain()
            query = {"DistributionOf": self._distribution_of.id_str()}

        elif values_of is not None:
            assert all(param is None for param in [data, known_terms, distribution_of])

            self._values_of = FindStatCompoundStatistic(values_of)
            domain = self._values_of.domain()
            query = {"ValuesOf": self._values_of.id_str()}

        else:
            raise ValueError("incompatible set of parameters: data: %s, distribution_of: %s, values_of: %s" % ((data, distribution_of, values_of)))

        if depth is not None:
            query["Depth"] = depth

        query["fields"] = "MatchingStatistic,Offset,Quality"
        if debug:
            print(query)
        verbose("querying FindStat %s" % query, caller_name='FindStatStatisticQuery')
        response = _post_json(FINDSTAT_API_STATISTICS, query)

        if debug:
            print(response)
        if "data" not in response:
            raise ValueError(response["error"])

        result = []
        for match in response["data"]:
            entry = response["included"]["MatchingStatistics"][match]
            result.append(FindStatMatchingStatistic(entry["MatchingStatistic"],
                                                    entry["Offset"],
                                                    entry["Quality"],
                                                    domain=domain))

        self._result = FancyTuple(result)

        FindStatFunction.__init__(self, FINDSTAT_STATISTIC_PADDED_IDENTIFIER % 0,
                                  data={"Bibliography": {},
                                        "Code": _get_code_from_callable(function),
                                        "Description": "",
                                        "Domain": domain,
                                        "Name": "a new statistic on %s" % domain.name("plural"),
                                        "References": "",
                                        "SageCode": ""},
                                  function=function)
        Element.__init__(self, FindStatStatistics()) # this is not completely correct, but it works

    def first_terms(self, max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        """
        Return the pairs of the known terms which contain singletons as a dictionary.

        EXAMPLES::

             sage: PM = PerfectMatchings
             sage: l = [(PM(2*n), [m.number_of_nestings() for m in PM(2*n)]) for n in range(5)]
             sage: r = findstat(l, depth=0); r                                  # optional -- internet
             0: St000041 (quality [99, 100])
             1: St000042 (quality [99, 100])
             sage: r.first_terms()                                              # optional -- internet
             {[]: 0, [(1, 2)]: 0}
        """
        return dict(itertools.islice(((objs[0], vals[0])
                                      for objs, vals in self._known_terms
                                      if len(vals) == 1), max_values))

    def _first_terms_raw(self, max_values):
        """
        Return the first terms as ``(string, value)`` pairs.

        EXAMPLES::

             sage: PM = PerfectMatchings
             sage: l = [(PM(2*n), [m.number_of_nestings() for m in PM(2*n)]) for n in range(5)]
             sage: r = findstat(l, depth=0); r                                  # optional -- internet
             0: St000041 (quality [99, 100])
             1: St000042 (quality [99, 100])
             sage: r._first_terms_raw(100)                                      # optional -- internet
             [('[]', 0), ('[(1, 2)]', 0)]
        """
        to_str = self.domain().to_string()
        return [(to_str(obj), val)
                for obj, val in self.first_terms(max_values=max_values).items()]

    def _generating_functions_dict(self,
                                   max_values=FINDSTAT_MAX_SUBMISSION_VALUES):
        """
        Return the generating functions of the levels where all values
        can be determined.

        TESTS::

            sage: n = 3; l = lambda i: [pi for pi in Permutations(n) if pi(1) == i]
            sage: data = [([pi for pi in l(i)], [pi(1) for pi in l(i)]) for i in range(1,n+1)]
            sage: data.append((Permutation([1,2]), 1))
            sage: q = findstat(data, depth=0); q                                # optional -- internet
            0: St000054 (quality [100, 100])
            sage: q.first_terms()                                               # optional -- internet
            {[1, 2]: 1}
            sage: q.generating_functions()                                      # optional -- internet, indirect doctest
            {3: 2*q^3 + 2*q^2 + 2*q}
        """
        return _distribution_from_data(self._known_terms,
                                       self.domain(),
                                       max_values,
                                       generating_functions=True)

    def __repr__(self):
        """
        Return a string representation of the query.

        EXAMPLES::

            sage: PM = PerfectMatchings
            sage: data = [(m, m.number_of_nestings()) for n in range(6) for m in PM(2*n)]
            sage: findstat(data, depth=1)                                       # optional -- internet
            0: St000042oMp00116 (quality [100, 100])
            1: St000041 (quality [20, 100])
            ...
        """
        if self._result:
            return repr(self._result)
        return "%s: %s" % (self.id_str(), self.name())

    def __getitem__(self, i):
        """
        Return the `t`-th result in the query.

        EXAMPLES::

            sage: PM = PerfectMatchings
            sage: data = [(m, m.number_of_nestings()) for n in range(6) for m in PM(2*n)]
            sage: r = findstat(data, depth=1)                                   # optional -- internet
            sage: r[1]                                                          # optional -- internet
            St000041 (quality [20, 100])
        """
        return self._result[i]

    def __len__(self):
        """
        Return the number of results in the query.

        EXAMPLES::

            sage: r = findstat(Permutations, lambda pi: pi.saliances()[0], depth=1)       # optional -- internet
            sage: len(r) > 4                                                              # optional -- internet
            True
        """
        return len(self._result)


class FindStatCompoundStatistic(Element, FindStatCombinatorialStatistic):
    def __init__(self, id, domain=None, check=True):
        """
        Initialize a compound statistic.

        A compound statistic is a sequence of maps followed by a statistic.

        INPUT:

        - ``id`` -- a padded identifier

        - ``domain`` -- (optional), the domain of the compound statistic

        - ``check`` -- whether to check that domains and codomains fit

        If the domain is given and ``check`` is ``False``, it is not
        fetched from FindStat.

        TESTS::

            sage: findstat("St000041oMp00127")                                  # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: the statistic St000041: The number of nestings of a perfect matching. cannot be composed with the map Mp00127
        """
        if isinstance(id, (int, Integer)):
            id = FINDSTAT_STATISTIC_PADDED_IDENTIFIER % id
        elif isinstance(id, FindStatCombinatorialStatistic):
            id = id.id_str()
        if domain is not None:
            self._domain = FindStatCollection(domain)
        else:
            self._domain = None
        composition = id.partition(FINDSTAT_MAP_SEPARATOR)
        self._statistic = FindStatStatistic(composition[0])
        if composition[2]:
            self._maps = FindStatCompoundMap(composition[2], domain=self._domain)
            self._id = self._statistic.id_str() + FINDSTAT_MAP_SEPARATOR + self._maps.id_str()
            if self._domain is None:
                self._domain = self._maps.domain()
        else:
            if self._domain is None:
                self._domain = self._statistic.domain()
            self._maps = FindStatCompoundMap("", domain=self._domain, codomain=self._domain)
            self._id = self._statistic.id_str()
        if (check
            and self._maps.codomain() != self._statistic.domain()):
            raise ValueError("the statistic %s cannot be composed with the map %s" % (self._statistic, self._maps))

        FindStatCombinatorialStatistic.__init__(self)
        Element.__init__(self, FindStatStatistics()) # this is not completely correct, but it works

    def _fetch_first_terms_raw(self):
        r"""
        Return the first terms of the compound statistic, as ``(string,
        value)`` pairs, fetched from FindStat.

        TESTS::

            sage: findstat("St000042oMp00116")._first_terms_raw(4)              # optional -- internet, indirect doctest
            [('[(1,2)]', 0),
             ('[(1,2),(3,4)]', 0),
             ('[(1,3),(2,4)]', 0),
             ('[(1,4),(2,3)]', 1)]
        """
        fields = "Values"
        url = FINDSTAT_API_STATISTICS + self.id_str() + "?fields=" + fields
        if len(self._maps):
            values = _get_json(url)["included"]["CompoundStatistics"][self.id_str()]["Values"]
        else:
            values = _get_json(url)["included"]["Statistics"][self.id_str()]["Values"]
        return [(sequence[0], sequence[-1]) for sequence in values]

    def domain(self):
        """
        Return the domain of the compound statistic.

        EXAMPLES::

            sage: findstat("St000042oMp00116").domain()                         # optional -- internet
            Cc0012: Perfect matchings
        """
        return self._domain

    def __call__(self, elt):
        """
        Apply the compound statistic to the given element.

        Note that this is only possible if execution of code is
        enabled, by setting the attribute ``_function`` of each map
        and the statistic to ``True``.

        EXAMPLES::

            sage: findstat("St000042oMp00116")(PerfectMatching([(1,2)]))        # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: execution of verified code provided by FindStat is not enabled for Mp00116: Kasraoui-Zeng
        """
        return self.statistic()(self.compound_map()(elt))

    def id_str(self):
        """
        Return the padded identifier of the compound statistic.

        EXAMPLES::

            sage: findstat("St000042oMp00116").id_str()                         # optional -- internet
            'St000042oMp00116'
        """
        return self._id

    def _repr_(self):
        """
        Return a string representation of the compound statistic.

        EXAMPLES::

            sage: findstat("St000042oMp00116")                                  # optional -- internet
            St000042oMp00116
        """
        return self.id_str()

    def statistic(self):
        """
        Return the statistic of the compound statistic.

        EXAMPLES::

            sage: findstat("St000041oMp00116").statistic()                      # optional -- internet
            St000041: The number of nestings of a perfect matching.
        """
        return self._statistic

    def compound_map(self):
        """
        Return the compound map which is part of the compound statistic.

        EXAMPLES::

            sage: findstat("St000051oMp00061oMp00069").compound_map()           # optional -- internet
            Mp00061oMp00069
        """
        return self._maps

    def browse(self):
        r"""
        Open the FindStat web page of the compound statistic in a browser.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCompoundStatistic
            sage: FindStatCompoundStatistic("St000042oMp00116").browse()        # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_STATISTICS + self.id_str())

    def info(self):
        """
        Print a detailed description of the compound statistic.

        EXAMPLES::

            sage: findstat("St000042oMp00116").info()                           # optional -- internet
                Mp00116: Kasraoui-Zeng: Perfect matchings -> Perfect matchings
                St000042: The number of crossings of a perfect matching.
        """
        if len(self.compound_map()):
            self.compound_map().info()
        self.statistic().info()


class FindStatMatchingStatistic(FindStatCompoundStatistic):
    def __init__(self, matching_statistic, offset, quality, domain=None):
        """
        Initialize a FindStat statistic match.

        INPUT:

        - ``matching_statistic`` -- a compound statistic identifier

        - ``offset`` -- the offset of the values, as provided by FindStat

        - ``quality`` -- the quality of the match, as provided by FindStat

        - ``domain`` -- (optional) the domain of the compound statistic

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingStatistic
            sage: FindStatMatchingStatistic("St000042oMp00116", 1, [17, 83])    # optional -- internet
            St000042oMp00116 with offset 1 (quality [17, 83])
        """
        self._quality = quality
        self._offset = offset
        # we can trust that matches have fitting domain / codomain sequence
        FindStatCompoundStatistic.__init__(self, matching_statistic, domain=domain, check=False)

    def _repr_(self):
        """
        Return a string representation of the match.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingStatistic
            sage: FindStatMatchingStatistic("St000042oMp00116", 1, [17, 83])    # optional -- internet
            St000042oMp00116 with offset 1 (quality [17, 83])
        """
        if self._offset:
            return "%s with offset %s (quality %s)" % (self.id_str(), self._offset, self._quality)
        return "%s (quality %s)" % (self.id_str(), self.quality())

    def offset(self):
        """
        Return the offset which has to be added to each value of the
        compound statistic to obtain the desired value.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingStatistic
            sage: r = FindStatMatchingStatistic("St000042oMp00116", 1, [17, 83])          # optional -- internet
            sage: r.offset()                                                    # optional -- internet
            1
        """
        return self._offset

    def quality(self):
        """
        Return the quality of the match, as provided by FindStat.

        The quality of a statistic match is a pair of percentages
        `(q_a, q_d)`, where `q_a` is the percentage of ``(object,
        value)`` pairs that are in the database among those which
        were sent to FindStat, and `q_d` is the percentage of
        ``(object, value)`` pairs with distinct values in the
        database among those which were sent to FindStat.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingStatistic
            sage: r = FindStatMatchingStatistic("St000042oMp00116", 1, [17, 83])          # optional -- internet
            sage: r.quality()                                                   # optional -- internet
            [17, 83]
        """
        return self._quality[:]

    def info(self):
        """
        Print a detailed explanation of the match.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingStatistic
            sage: r = FindStatMatchingStatistic("St000042oMp00116", 1, [17, 83])          # optional -- internet
            sage: r.info()                                                      # optional -- internet
            after adding 1 to every value
            and applying
                Mp00116: Kasraoui-Zeng: Perfect matchings -> Perfect matchings
            to the objects (see `.compound_map()` for details)
            <BLANKLINE>
            your input matches
                St000042: The number of crossings of a perfect matching.
            <BLANKLINE>
            among the values you sent, 17 percent are actually in the database,
            among the distinct values you sent, 83 percent are actually in the database

            sage: r = FindStatMatchingStatistic("St000042", 1, [17, 83])        # optional -- internet
            sage: r.info()                                                      # optional -- internet
            after adding 1 to every value
            <BLANKLINE>
            your input matches
                St000042: The number of crossings of a perfect matching.
            <BLANKLINE>
            among the values you sent, 17 percent are actually in the database,
            among the distinct values you sent, 83 percent are actually in the database
        """
        if self.offset() < 0:
            print("after subtracting %s from every value" % (-self.offset()))
        if self.offset() > 0:
            print("after adding %s to every value" % self.offset())
        if len(self.compound_map()):
            if self.offset():
                print("and applying")
            else:
                print("after applying")
            self.compound_map().info()
            print("to the objects (see `.compound_map()` for details)")
        print()
        print("your input matches")
        self.statistic().info()
        print()
        print("among the values you sent, %s percent are actually in the database," % self.quality()[0])
        print("among the distinct values you sent, %s percent are actually in the database" % self.quality()[1])

######################################################################
# maps
######################################################################


class FindStatCombinatorialMap(SageObject):
    """
    A class serving as common ancestor of :class:`FindStatStatistic`
    and :class:`FindStatCompoundStatistic`.
    """
    pass


class FindStatMap(Element,
                  FindStatFunction,
                  FindStatCombinatorialMap,
                  metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A FindStat map.

    :class:`FindStatMap` is a class representing a combinatorial
    map available in the FindStat database.

    This class provides methods to inspect various properties of
    these maps, in particular :meth:`code`.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatMap
        sage: FindStatMap(116)                                                  # optional -- internet
        Mp00116: Kasraoui-Zeng

    .. SEEALSO::

        :class:`FindStatMaps`
    """
    @staticmethod
    def __classcall_private__(cls, entry):
        """
        Retrieve a map from the database.

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap("abcdefgh")                                       # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: the value 'abcdefgh' is not a valid FindStat map identifier
        """
        return FindStatMaps()(entry)

    def __init__(self, parent, id):
        """
        Initialize the map.

        This should only be called in
        :meth:`FindStatMaps()._element_constructor_` via
        :meth:`FindStatMaps().element_class`.

        INPUT:

        - ``parent`` -- :class:`FindStatMaps`

        - ``id`` -- the (padded) FindStat identifier of the statistic

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(116).parent()                                     # optional -- internet
            Set of combinatorial maps used by FindStat
        """
        FindStatFunction.__init__(self, id)
        Element.__init__(self, parent)

    def __reduce__(self):
        """
        Return a function and its arguments needed to create the map.

        TESTS::

            sage: from sage.databases.findstat import FindStatMap
            sage: c = FindStatMap(116)                                          # optional -- internet
            sage: loads(dumps(c)) == c                                          # optional -- internet
            True
        """
        return (FindStatMap, (self.id(),))

    def _richcmp_(self, other, op):
        """
        Compare two maps by identifier.

        TESTS::

            sage: findmap(61) != findmap(62)                                    # optional -- internet
            True
            sage: findmap(61) == findstat(61)                                   # optional -- internet
            False
        """
        return richcmp(self.id(), other.id(), op)

    def _fetch_data(self):
        r"""
        Return a dictionary containing the data of the map, fetched from
        FindStat.

        TESTS::

            sage: findmap(64)._data()                                           # optional -- internet, indirect doctest
            {'Bibliography': {},
             'Codomain': 'Cc0001',
             'Description': 'Sends a permutation to its reverse.\r\n\r\nThe reverse of a permutation $\\sigma$ of length $n$ is given by $\\tau$ with $\\tau(i) = \\sigma(n+1-i)$.',
             'Domain': 'Cc0001',
             'Name': 'reverse',
             'Properties': 'bijective, graded, involutive',
             'References': '',
             'SageCode': 'def mapping(sigma):\r\n    return sigma.reverse()'}
        """
        fields = "Bibliography,Codomain,Description,Domain,Name,Properties,References,SageCode"
        fields_Bibliography = "Author,Title"
        url = (FINDSTAT_API_MAPS + self.id_str()
               + "?fields=" + fields
               + "&fields[Bibliography]=" + fields_Bibliography)
        verbose("fetching map data %s" % url, caller_name='FindStatMap')
        included = _get_json(url)["included"]
        # slightly simplify the representation
        data = included["Maps"][self.id_str()]
        # we replace the list of identifiers in Bibliography with the dictionary
        data["Bibliography"] = included["References"]
        return data

    def browse(self):
        r"""
        Open the FindStat web page of the map in a browser.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(116).browse()                                     # optional -- webbrowser
        """
        if self.id() == 0:
            self.submit()
        else:
            webbrowser.open(FINDSTAT_URL_MAPS + self.id_str())

    def submit(self):
        r"""
        Open the FindStat web page for editing the map in a browser.

        TESTS::

            sage: s = findmap(62)                                               # optional -- internet
            sage: s.set_name(u"Möbius")                                         # optional -- internet
            sage: s.submit()                                                    # optional -- webbrowser
            sage: s.reset()                                                     # optional -- internet
        """
        args = dict()
        args["OriginalMap"] = self.id_str()
        args["Domain"] = self.domain().id_str()
        args["Codomain"] = self.codomain().id_str()
        args["Name"] = self.name()
        args["Description"] = self.description()
        args["References"] = self.references_raw()
        args["Properties"] = self.properties_raw()
        args["SageCode"] = self.sage_code()
        args["CurrentAuthor"] = FindStat().user_name()
        args["CurrentEmail"] = FindStat().user_email()

        if not self.id():
            url = FINDSTAT_NEWMAP_FORM_HEADER % FINDSTAT_URL_NEW_MAP
        else:
            url = FINDSTAT_NEWMAP_FORM_HEADER % (FINDSTAT_URL_EDIT_MAP + self.id_str())
        _submit(args, url)

    # editing and submitting is really the same thing
    edit = submit

    def __hash__(self):
        """
        Return a hash value for the map.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: sorted(list(FindStatMaps(domain=1, codomain=10)))             # optional -- internet, indirect doctest
            [Mp00061: to increasing tree, Mp00072: binary search tree: left to right]
        """
        return self.id()

    def codomain(self):
        r"""
        Return the FindStat collection which is the codomain of the map.

        OUTPUT: the codomain of the map as a :class:`FindStatCollection`

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: FindStatMap(27).codomain()                                    # optional -- internet
            Cc0002: Integer partitions
        """
        return FindStatCollection(self._data()["Codomain"])

    def properties_raw(self):
        r"""
        Return the properties of the map.

        OUTPUT: the properties as a string

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap               # optional -- internet
            sage: FindStatMap(61).properties_raw()                              # optional -- internet
            'surjective, graded'
        """
        return self._data()["Properties"]

    def set_properties_raw(self, value):
        r"""
        Set the properties of the map.

        EXAMPLES::

            sage: # optional - internet
            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(61).set_properties_raw('surjective')
            sage: FindStatMap(61).properties_raw()
            'surjective'
            sage: FindStatMap(61)
            Mp00061(modified): to increasing tree
            sage: FindStatMap(61).reset()
            sage: FindStatMap(61)
            Mp00061: to increasing tree
        """
        if value != self.properties_raw():
            self._modified = True
            self._data_cache["Properties"] = value

    def set_name(self, value):
        r"""
        Set the name of the map.

        INPUT:

        - ``value`` -- string; the new name of the map

        This information is used when submitting the map with
        :meth:`submit`.

        TESTS::

            sage: s = findmap(62)                                               # optional -- internet
            sage: s.set_name(u"Möbius"); s                                      # optional -- internet
            Mp00062(modified): Möbius
            sage: s.reset(); s                                                  # optional -- internet
            Mp00062: Lehmer-code to major-code bijection
        """
        if value != self.name():
            self._modified = True
            self._data_cache["Name"] = value

    def info(self):
        """
        Print a detailed description of the map.

        EXAMPLES::

            sage: findmap("Mp00116").info()                                     # optional -- internet
                Mp00116: Kasraoui-Zeng: Perfect matchings -> Perfect matchings
        """
        print("    %s: %s -> %s" % (self,
                                    self.domain().name("plural"),
                                    self.codomain().name("plural")))


_all_maps = {}


class FindStatMaps(UniqueRepresentation, Parent):
    r"""
    The class of FindStat maps.

    The elements of this class are combinatorial maps currently in
    FindStat.

    EXAMPLES:

    We can print a sample map for each domain and codomain::

        sage: from sage.databases.findstat import FindStatCollections, FindStatMaps
        sage: ccs = sorted(FindStatCollections())[:3]                           # optional -- internet
        sage: for cc_dom in ccs:                                                # optional -- internet
        ....:     for cc_codom in ccs:
        ....:         print(cc_dom.name(style='plural') + " -> " + cc_codom.name(style='plural'))
        ....:         try:
        ....:             print("    " + next(iter(FindStatMaps(cc_dom, cc_codom))).name())
        ....:         except StopIteration:
        ....:             pass
        Permutations -> Permutations
            Lehmer-code to major-code bijection
        Permutations -> Integer partitions
            Robinson-Schensted tableau shape
        Permutations -> Dyck paths
            left-to-right-maxima to Dyck path
        Integer partitions -> Permutations
        Integer partitions -> Integer partitions
            conjugate
        Integer partitions -> Dyck paths
            to Dyck path
        Dyck paths -> Permutations
            to non-crossing permutation
        Dyck paths -> Integer partitions
            to partition
        Dyck paths -> Dyck paths
            reverse
    """
    def __init__(self, domain=None, codomain=None):
        """
        TESTS::

            sage: from sage.databases.findstat import FindStatMaps
            sage: M = FindStatMaps()                                            # optional -- internet
            sage: TestSuite(M).run()                                            # optional -- internet
        """
        if domain is None:
            self._domain = None
        else:
            self._domain = FindStatCollection(domain)
        if codomain is None:
            self._codomain = None
        else:
            self._codomain = FindStatCollection(codomain)
        self._identifiers = None
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, id):
        """
        Initialize a FindStat map.

        INPUT:

        - ``id`` -- string containing the FindStat identifier of
          the map, or an integer giving its id

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMap
            sage: FindStatMap(61)                                               # optional -- internet, indirect doctest
            Mp00061: to increasing tree
        """
        if isinstance(id, self.Element):
            return id
        if isinstance(id, (int, Integer)):
            id = FINDSTAT_MAP_PADDED_IDENTIFIER % id
        elif isinstance(id, FindStatCombinatorialMap):
            id = id.id_str()
        if not isinstance(id, str):
            raise TypeError("the value %s is neither an integer nor a string" % id)
        else:
            id = id.strip()
        if FINDSTAT_MAP_SEPARATOR in id:
            return FindStatCompoundMap(id)
        if not re.match(FINDSTAT_MAP_REGEXP, id) or id == FINDSTAT_MAP_PADDED_IDENTIFIER % 0:
            raise ValueError("the value '%s' is not a valid FindStat map identifier" % id)
        if id not in _all_maps or _all_maps[id] is None:
            _all_maps[id] = self.element_class(self, id)

        return _all_maps[id]

    def _repr_(self):
        """
        Return the representation of the set of FindStat maps.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: FindStatMaps()                                                # optional -- internet
            Set of combinatorial maps used by FindStat
        """
        if self._domain is None:
            text = []
        else:
            text = ["domain %s" % self._domain]
        if self._codomain is not None:
            text.append("codomain %s" % self._codomain)

        if text:
            return "Set of combinatorial maps with " + " and ".join(text) + " used by FindStat"
        return "Set of combinatorial maps used by FindStat"

    def __iter__(self):
        """
        Return an iterator over all FindStat maps.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMaps
            sage: next(iter(FindStatMaps(domain=1, codomain=10)))               # optional -- internet
            Mp00061: to increasing tree
        """
        if self._identifiers is None:
            query = []
            if self._domain is not None:
                query.append("Domain=%s" % self._domain.id_str())
            if self._codomain is not None:
                query.append("Codomain=%s" % self._codomain.id_str())
            if query:
                url = FINDSTAT_API_MAPS + "?" + "&".join(query)
            else:
                url = FINDSTAT_API_MAPS
            self._identifiers = _get_json(url)["data"]

        for mp in self._identifiers:
            yield FindStatMap(mp)

    def _an_element_(self):
        """
        Return a FindStat map.

        EXAMPLES::

            sage: findmap(domain="Dyck paths", codomain='Posets').an_element()  # optional -- internet
            Mp00232: parallelogram poset
        """
        try:
            return next(iter(self))
        except StopIteration:
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError

    Element = FindStatMap


class FindStatMapQuery(FindStatMap):
    """
    A class representing a query for FindStat (compound) maps.
    """
    def __init__(self, data=None, values_of=None, distribution_of=None,
                 domain=None, codomain=None, known_terms=None, function=None,
                 depth=FINDSTAT_DEFAULT_DEPTH,
                 debug=False):
        """
        Initialize a query for FindStat (compound) maps.

        INPUT:

        - ``data`` -- (optional), a list of pairs ``(objects,
          values)``, where ``objects`` and ``values`` are all lists
          of the same length, the former are elements in the domain
          and the latter in the codomain

        - ``known_terms`` -- (optional), a lazy list in the same format
          as ``data``, which agrees with ``data``, and may be used
          for submission

        - ``values_of`` -- (optional), anything accepted by
          :class:`FindStatCompoundStatistic`

        - ``distribution_of`` -- (optional), anything accepted by
          :class:`FindStatCompoundStatistic`

        - ``domain``, ``codomain`` -- (optional), anything accepted by
          :class:`FindStatCollection`

        - ``depth`` -- (optional), the number of maps to apply before
          applying the statistic

        - ``function`` -- (optional), a callable producing the terms

        Only one of ``data``, ``values_of`` and ``distribution_of``
        may be provided.  The parameter ``domain`` must be provided
        if and only if ``data`` is provided, or ``values_of`` or
        ``distribution_of`` are given as a function.

        The parameter ``known_terms`` is only allowed, if ``data`` is
        provided.  It defaults to ``data``.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMapQuery
            sage: data = [[[pi], [pi.complement().increasing_tree_shape()]] for pi in Permutations(4)]
            sage: FindStatMapQuery(domain=1, codomain=10, data=data)            # optional -- internet
            0: Mp00061oMp00069 (quality [100])
        """
        self._first_terms = data
        if data is not None and known_terms is None:
            self._known_terms = data
        else:
            self._known_terms = known_terms
        self._values_of = None
        self._distribution_of = None
        self._depth = depth

        if data is not None:
            assert all(param is None for param in [distribution_of, values_of])

            domain = FindStatCollection(domain)
            codomain = FindStatCollection(codomain)
            query = {"Domain": domain.id_str(),
                     "Codomain": codomain.id_str(),
                     "Data": _data_to_str(self._first_terms, domain, codomain)}

        elif distribution_of is not None:
            assert all(param is None for param in [data, known_terms, values_of])

            self._distribution_of = FindStatCompoundMap(distribution_of)
            domain = self._distribution_of.domain()
            codomain = self._distribution_of.codomain()
            query = {"DistributionOf": self._distribution_of.id_str()}

        elif values_of is not None:
            assert all(param is None for param in [data, known_terms, distribution_of])

            self._values_of = FindStatCompoundMap(values_of)
            domain = self._values_of.domain()
            codomain = self._values_of.codomain()
            query = {"ValuesOf": self._values_of.id_str()}

        else:
            raise ValueError("incompatible set of parameters: data: %s, distribution_of: %s, values_of: %s" % ((data, distribution_of, values_of)))

        if depth is not None:
            query["Depth"] = depth

        query["fields"] = "MatchingMap,Quality"
        if debug:
            print(query)
        verbose("querying FindStat %s" % query, caller_name='FindStatMapQuery')
        response = _post_json(FINDSTAT_API_MAPS, query)

        if debug:
            print(response)
        if "data" not in response:
            raise ValueError(response["error"])

        result = []
        for match in response["data"]:
            entry = response["included"]["MatchingMaps"][match]
            result.append(FindStatMatchingMap(entry["MatchingMap"],
                                              entry["Quality"]))
        self._result = FancyTuple(result)

        FindStatFunction.__init__(self, FINDSTAT_MAP_PADDED_IDENTIFIER % 0,
                                  data={"Bibliography": {},
                                        "Code": _get_code_from_callable(function),
                                        "Description": "",
                                        "Domain": domain,
                                        "Codomain": codomain,
                                        "Name": "a new map from %s to %s" % (domain.name("plural"), codomain.name("plural")),
                                        "References": "",
                                        "Properties": "",
                                        "SageCode": ""},
                                  function=function)
        Element.__init__(self, FindStatMaps()) # this is not completely correct, but it works

    def __repr__(self):
        """
        Return a string representation of the query.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMapQuery
            sage: data = [[[pi],[pi.complement().increasing_tree_shape()]] for pi in Permutations(4)]
            sage: FindStatMapQuery(domain=1, codomain=10, data=data)            # optional -- internet
            0: Mp00061oMp00069 (quality [100])
        """
        if self._result:
            return repr(self._result)
        return "%s: %s" % (self.id_str(), self.name())

    def __getitem__(self, i):
        """
        Return the `i`-th result in the query.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMapQuery
            sage: data = [[[pi],[pi.complement().increasing_tree_shape()]] for pi in Permutations(4)]
            sage: r = FindStatMapQuery(domain=1, codomain=10, data=data)        # optional -- internet
            sage: r[0]                                                          # optional -- internet
            Mp00061oMp00069 (quality [100])
        """
        return self._result[i]

    def __len__(self):
        """
        Return the number of results in the query.

        EXAMPLES::

            sage: r = findmap("Permutations", lambda pi: pi.descents_composition())       # optional -- internet
            sage: len(r) > 0                                                              # optional -- internet
            True
        """
        return len(self._result)


class FindStatCompoundMap(Element, FindStatCombinatorialMap):
    def __init__(self, id, domain=None, codomain=None, check=True):
        """
        Initialize a compound statistic.

        INPUT:

        - ``id`` -- a padded identifier

        - ``domain`` -- (optional), the domain of the compound map

        - ``codomain`` -- (optional), the codomain of the compound map

        - ``check`` -- whether to check that domains and codomains fit

        If domain and codomain are given and ``check`` is ``False``,
        they are not fetched from FindStat.

        If ``id`` is the empty string, ``domain`` must be provided,
        and the identity map on this collection is returned.

        TESTS::

            sage: findmap("Mp00146oMp00127").domain()                           # optional -- internet
            Cc0001: Permutations

            sage: findmap("Mp00146oMp00127").codomain()                         # optional -- internet
            Cc0012: Perfect matchings

            sage: findmap("Mp00127oMp00146")                                    # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: the sequence of maps [Mp00146: to tunnel matching, Mp00127: left-to-right-maxima to Dyck path] cannot be composed
        """
        if isinstance(id, (int, Integer)):
            id = FINDSTAT_MAP_PADDED_IDENTIFIER % id
        elif isinstance(id, FindStatCombinatorialMap):
            id = id.id_str()
        if id == "":
            self._id = "id"
            self._maps = []
            self._domain = FindStatCollection(domain)
            self._codomain = self._domain
        else:
            self._maps = [FindStatMap(m) for m in id.split(FINDSTAT_MAP_SEPARATOR)][::-1]
            if (check
                and not all(self._maps[i].codomain() == self._maps[i+1].domain()
                            for i in range(len(self._maps)-1))):
                raise ValueError("the sequence of maps %s cannot be composed" % self._maps)
            if domain is None:
                self._domain = self._maps[0].domain()
            else:
                self._domain = FindStatCollection(domain)
            if codomain is None:
                self._codomain = self._maps[-1].codomain()
            else:
                self._codomain = FindStatCollection(codomain)
            self._id = FINDSTAT_MAP_SEPARATOR.join(m.id_str() for m in reversed(self._maps))

        Element.__init__(self, FindStatMaps()) # this is not completely correct, but it works

    def domain(self):
        """
        Return the domain of the compound map.

        EXAMPLES::

            sage: findmap("Mp00099oMp00127").domain()                           # optional -- internet
            Cc0001: Permutations
        """
        return self._domain

    def codomain(self):
        """
        Return the codomain of the compound map.

        EXAMPLES::

            sage: findmap("Mp00099oMp00127").codomain()                         # optional -- internet
            Cc0005: Dyck paths
        """
        return self._codomain

    def __call__(self, elt):
        """
        Apply the compound map to the given element.

        Note that this is only possible if execution of code is
        enabled, by setting the attribute ``_function`` of each map
        to ``True``.

        EXAMPLES::

            sage: findmap("Mp00099oMp00127")(Permutation([1,3,2]))              # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: execution of verified code provided by FindStat is not enabled for Mp00127: left-to-right-maxima to Dyck path
        """
        for m in self.maps():
            elt = m(elt)
        return elt

    def __getitem__(self, i):
        """
        Return the `i`-th map in the compound map.

        EXAMPLES::

            sage: findmap("Mp00099oMp00127")[1]                                 # optional -- internet
            Mp00099: bounce path
        """
        return self._maps[i]

    def id_str(self):
        """
        Return the padded identifier of the compound map.

        EXAMPLES::

            sage: findmap("Mp00099oMp00127").id_str()                           # optional -- internet
            'Mp00099oMp00127'
        """
        return self._id

    def _repr_(self):
        """
        Return a string representation of the compound statistic.

        EXAMPLES::

            sage: findmap("Mp00099oMp00127")                                    # optional -- internet
            Mp00099oMp00127
        """
        return self.id_str()

    def maps(self):
        """
        Return the maps occurring in the compound map as a list.

        EXAMPLES::

            sage: findmap("Mp00099oMp00127").maps()                             # optional -- internet
            [Mp00127: left-to-right-maxima to Dyck path, Mp00099: bounce path]
        """
        return self._maps

    def __len__(self):
        """
        Return the number of maps occurring in the compound map.

        .. WARNING::

            Do not confuse this with the number of results in a query.

        EXAMPLES::

            sage: len(findmap("Mp00099oMp00127"))                               # optional -- internet
            2
        """
        return len(self._maps)

    def browse(self):
        r"""
        Open the FindStat web page of the compound map in a browser.

        EXAMPLES::

            sage: findmap(62).browse()                                          # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_MAPS + self.id_str())

    def info(self):
        """
        Print a detailed explanation of the compound map.

        EXAMPLES::

            sage: findmap("Mp00099oMp00127").info()                             # optional -- internet
                Mp00127: left-to-right-maxima to Dyck path: Permutations -> Dyck paths
                Mp00099: bounce path: Dyck paths -> Dyck paths
        """
        for mp in self:
                mp.info()


class FindStatMatchingMap(FindStatCompoundMap):
    def __init__(self, matching_map, quality, domain=None, codomain=None):
        """
        Initialize a FindStat map match.

        INPUT:

        - ``matching_map`` -- a compound map identifier

        - ``quality`` -- the quality of the match, as provided by FindStat

        - ``domain`` -- (optional) the domain of the compound map

        - ``codomain`` -- (optional), the codomain of the compound map

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingMap
            sage: FindStatMatchingMap("Mp00099oMp00127", [83])                  # optional -- internet
            Mp00099oMp00127 (quality [83])
        """
        self._quality = quality
        # we can trust that matches have fitting domain / codomain sequence
        FindStatCompoundMap.__init__(self, matching_map, domain=domain, codomain=codomain, check=False)

    def _repr_(self):
        """
        Return a string representation of the match.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingMap
            sage: FindStatMatchingMap("Mp00099oMp00127", [83])                  # optional -- internet
            Mp00099oMp00127 (quality [83])
        """
        return "%s (quality %s)" % (self.id_str(), self.quality())

    def quality(self):
        """
        Return the quality of the match, as provided by FindStat.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingMap
            sage: FindStatMatchingMap("Mp00099oMp00127", [83]).quality()        # optional -- internet
            [83]
        """
        return self._quality[:]

    def info(self):
        """
        Print a detailed explanation of the match.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatMatchingMap
            sage: FindStatMatchingMap("Mp00099oMp00127", [83]).info()           # optional -- internet
            your input matches
                Mp00127: left-to-right-maxima to Dyck path: Permutations -> Dyck paths
                Mp00099: bounce path: Dyck paths -> Dyck paths
            <BLANKLINE>
            among the values you sent, 83 percent are actually in the database
        """
        print("your input matches")
        super().info()
        print()
        print("among the values you sent, %s percent are actually in the database" % self.quality()[0])


######################################################################
# collections
######################################################################

# helper for generation of CartanTypes
def _finite_irreducible_cartan_types_by_rank(n):
    """
    Return the Cartan types of rank `n`.

    INPUT:

    - ``n`` -- integer

    OUTPUT: the list of Cartan types of rank `n`

    TESTS::

        sage: from sage.databases.findstat import _finite_irreducible_cartan_types_by_rank
        sage: _finite_irreducible_cartan_types_by_rank(2)
        [['A', 2], ['B', 2], ['G', 2]]
    """
    cartan_types = [CartanType(['A',n])]
    if n >= 2:
        cartan_types += [CartanType(['B',n])]
    if n >= 3:
        cartan_types += [CartanType(['C',n])]
    if n >= 4:
        cartan_types += [CartanType(['D',n])]
    if 6 <= n <= 8:
        cartan_types += [CartanType(['E',n])]
    if n == 4:
        cartan_types += [CartanType(['F',n])]
    if n == 2:
        cartan_types += [CartanType(['G',n])]
    return cartan_types

# helper for generation of PlanePartitions


def _plane_partitions_by_size_aux(n, outer=None):
    """
    Iterate over the plane partitions with `n` boxes, as lists.

    INPUT:

    - ``n`` -- integer

    OUTPUT: the plane partitions with `n` boxes as lists

    TESTS::

        sage: from sage.databases.findstat import _plane_partitions_by_size_aux
        sage: list(_plane_partitions_by_size_aux(3))
        [[[1], [1], [1]], [[2], [1]], [[1, 1], [1]], [[3]], [[2, 1]], [[1, 1, 1]]]
    """
    if n == 0:
        yield []
        return
    if outer is None:
        outer = [n]*n
    for k in range(1, n+1):
        for la in Partitions(k, outer=outer):
            for pp in _plane_partitions_by_size_aux(n-k, outer=la):
                pp = [la] + pp
                yield pp


def _plane_partitions_by_size(n):
    """
    Iterate over the plane partitions with `n` boxes.

    .. TODO::

        This can be replaced when :issue:`28244` is merged.

    INPUT:

    - ``n`` -- integer

    OUTPUT: the plane partitions with `n` boxes

    TESTS::

        sage: from sage.databases.findstat import _plane_partitions_by_size
        sage: list(_plane_partitions_by_size(3))
        [Plane partition [[1], [1], [1]],
         Plane partition [[2], [1]],
         Plane partition [[1, 1], [1]],
         Plane partition [[3]],
         Plane partition [[2, 1]],
         Plane partition [[1, 1, 1]]]
    """
    for pp in _plane_partitions_by_size_aux(n):
        yield PlanePartition(pp)

# helper for generation of Lattices


def _finite_lattices(n):
    """
    Iterate over the lattices with `n` elements.

    INPUT:

    - ``n`` -- integer

    OUTPUT: the lattices with `n` elements

    TESTS::

        sage: from sage.databases.findstat import _finite_lattices
        sage: [L.cover_relations() for L in _finite_lattices(4)]
        [[['bottom', 0], ['bottom', 1], [0, 'top'], [1, 'top']],
         [['bottom', 0], [0, 1], [1, 'top']]]
    """
    if n <= 2:
        for P in Posets(n):
            if P.is_lattice():
                yield LatticePoset(P)
    else:
        for P in Posets(n-2):
            Q = P.with_bounds()
            if Q.is_lattice():
                yield LatticePoset(Q)


class FindStatCollection(Element,
                         metaclass=InheritComparisonClasscallMetaclass):
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

    - a Sage object belonging to a collection, or

    - an iterable producing a Sage object belonging to a collection.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatCollection
        sage: FindStatCollection("Dyck paths")                                  # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(5)                                             # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(DyckWord([1,0,1,0]))                           # optional -- internet
        Cc0005: Dyck paths

        sage: FindStatCollection(DyckWords(2))                                  # optional -- internet
        a subset of Cc0005: Dyck paths

        sage: FindStatCollection(DyckWords)                                     # optional -- internet
        Cc0005: Dyck paths

    .. SEEALSO::

        :class:`FindStatCollections`
    """
    @staticmethod
    def __classcall_private__(cls, entry):
        """
        Retrieve a collection from the database.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection(0)                                         # optional -- internet
            Traceback (most recent call last):
            ...
            ValueError: could not find FindStat collection for 0
        """
        return FindStatCollections()(entry)

    def __init__(self, parent, id, data, sageconstructor_overridden):
        """
        Initialize the collection.

        This should only be called in
        :meth:`FindStatCollections()._element_constructor_` via
        :meth:`FindStatCollections().element_class`.

        INPUT:

        - ``parent`` -- :class:`FindStatCollections`

        - ``id`` -- the (padded) FindStat identifier of the collection

        - ``data`` -- dictionary containing the properties of the
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
        self._data = data
        self._sageconstructor_overridden = sageconstructor_overridden

        Element.__init__(self, parent)

    def __reduce__(self):
        """
        Return a function and its arguments needed to create the
        collection.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("Permutations")                        # optional -- internet
            sage: loads(dumps(c)) == c                                          # optional -- internet
            True
        """
        return FindStatCollection, (self.id(),)

    def _richcmp_(self, other, op):
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

        It is not clear what we want here, but equality may be better::

            sage: FindStatCollection(graphs(3)) == FindStatCollection("Graphs")                     # optional -- internet
            True
        """
        return richcmp(self.id(), other.id(), op)

    def __hash__(self):
        """
        Return a hash value for the collection.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollections
            sage: set(FindStatCollections())           # optional -- internet
            {Cc0001: Permutations,
             Cc0002: Integer partitions,
             ...
        """
        return self.id()

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
        return "Code" in self._data

    def elements_on_level(self, level):
        """
        Return an iterable over the elements on the given level.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Perfect Matchings").elements_on_level(4)  # optional -- internet
            Perfect matchings of {1, 2, 3, 4}
        """
        return self._data["Code"].elements_on_level(level)

    def element_level(self, element):
        """
        Return the level of an element.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: cc = FindStatCollection("Perfect Matchings")                  # optional -- internet
            sage: cc.element_level(PerfectMatching([[1,2],[3,4],[5,6]]))        # optional -- internet
            6
        """
        return self._data["Code"].element_level(element)

    def is_element(self, element):
        """
        Return whether the element belongs to the collection.

        If the collection is not yet supported, return whether element is a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: cc = FindStatCollection("Perfect Matchings")                  # optional -- internet
            sage: cc.is_element(PerfectMatching([[1,2],[3,4],[5,6]]))           # optional -- internet
            True

            sage: cc.is_element(SetPartition([[1,2],[3,4],[5,6]]))              # optional -- internet
            False
        """
        if self.is_supported():
            return self._data["Code"].is_element(element)

        return isinstance(element, str)

    def levels_with_sizes(self):
        """
        Return a dictionary from levels to level sizes.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: cc = FindStatCollection("Perfect Matchings")                  # optional -- internet
            sage: cc.levels_with_sizes()                                        # optional -- internet
            {2: 1, 4: 3, 6: 15, 8: 105, 10: 945}
        """
        return self._data["LevelsWithSizes"]

    def in_range(self, element):
        r"""
        Check whether an element of the collection is in FindStat's precomputed range.

        INPUT:

        - ``element`` -- a Sage object that belongs to the collection

        OUTPUT: ``True``, if ``element`` is used by the FindStat search
        engine, and ``False`` if it is ignored

        EXAMPLES::

            sage: # optional - internet
            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")
            sage: c.in_range(GelfandTsetlinPattern([[2, 1], [1]]))
            True
            sage: c.in_range(GelfandTsetlinPattern([[3, 1], [1]]))
            True
            sage: c.in_range(GelfandTsetlinPattern([[7, 1], [1]]))
            False

        TESTS::

            sage: from sage.databases.findstat import FindStatCollections
            sage: l = FindStatCollections()                                     # optional -- internet
            sage: long = [14, 20]
            sage: for c in sorted(l):                                           # optional -- internet
            ....:     if c.id() not in long and c.is_supported():
            ....:         f = c.first_terms(lambda x: 1)
            ....:         print("{} {} {}".format(c, len(list(f)), all(c.in_range(e) for e, _ in f)))
            Cc0001: Permutations 5913 True
            Cc0002: Integer partitions 1211 True
            Cc0005: Dyck paths 2055 True
            Cc0006: Integer compositions 1023 True
            Cc0007: Standard tableaux 1115 True
            Cc0009: Set partitions 1155 True
            Cc0010: Binary trees 2055 True
            Cc0012: Perfect matchings 1069 True
            Cc0013: Cores 100 True
            Cc0017: Alternating sign matrices 7917 True
            Cc0018: Gelfand-Tsetlin patterns 1409 True
            Cc0019: Semistandard tableaux 2374 True
            Cc0021: Ordered trees 2056 True
            Cc0022: Finite Cartan types 31 True
            Cc0023: Parking functions 18248 True
            Cc0024: Binary words 1022 True
            Cc0025: Plane partitions 1123 True
            Cc0026: Decorated permutations 2371 True
            Cc0027: Signed permutations 4282 True
            Cc0028: Skew partitions 1250 True
            Cc0029: Lattices 1378 True
            Cc0030: Ordered set partitions 5316 True
        """
        return self._data["Code"].element_level(element) in self._data["LevelsWithSizes"]

    def first_terms(self, function, level=None):
        r"""
        Compute the first few terms of the given function, possibly
        restricted to a level, as a lazy list.

        INPUT:

        - ``function`` -- a callable

        - ``level`` -- (optional), the level to restrict to

        OUTPUT:

        A lazy list of pairs of the form ``(object, value)``.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.first_terms(lambda x: 1)[:10].list()                        # optional -- internet
            [([[1, 0], [0]], 1),
             ([[1, 0], [1]], 1),
             ([[2, 0], [0]], 1),
             ([[2, 0], [1]], 1),
             ([[2, 0], [2]], 1),
             ([[1, 1], [1]], 1),
             ([[1, 0, 0], [0, 0], [0]], 1),
             ([[1, 0, 0], [1, 0], [0]], 1),
             ([[1, 0, 0], [1, 0], [1]], 1),
             ([[3, 0], [0]], 1)]
        """
        if self._sageconstructor_overridden is None:
            if level is None:
                g = (x
                     for level in self._data["LevelsWithSizes"]
                     for x in self._data["Code"].elements_on_level(level))
            else:
                g = (x for x in self._data["Code"].elements_on_level(level))
        else:
            if level is None:
                g = self._sageconstructor_overridden
            else:
                g = (x for x in self._sageconstructor_overridden
                     if self.element_level(x) == level)

        return lazy_list(((x, function(x)) for x in g))

    def id(self):
        r"""
        Return the FindStat identifier of the collection.

        OUTPUT: the FindStat identifier of the collection as an integer

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.id()                                                        # optional -- internet
            18
        """
        return int(self._id[2:])

    def id_str(self):
        r"""
        Return the FindStat identifier of the collection.

        OUTPUT: the FindStat identifier of the collection as a string

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("GelfandTsetlinPatterns")              # optional -- internet
            sage: c.id_str()                                                    # optional -- internet
            'Cc0018'
        """
        return self._id

    def browse(self):
        r"""
        Open the FindStat web page of the collection in a browser.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Permutations").browse()                   # optional -- webbrowser
        """
        webbrowser.open(FINDSTAT_URL_COLLECTIONS + self._id)

    def to_string(self):
        r"""
        Return a function that returns a FindStat representation given an
        object.

        If the collection is not yet supported, return the identity.

        OUTPUT:

        The function that produces the string representation as
        needed by the FindStat search webpage.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: p = Poset((range(3), [[0, 1], [1, 2]]))                       # optional -- internet
            sage: c = FindStatCollection("Posets")                              # optional -- internet
            sage: c.to_string()(p)                                              # optional -- internet
            '([(0, 1), (1, 2)], 3)'
        """
        if self.is_supported():
            return self._data["Code"].element_to_string
        return lambda x: x

    def from_string(self):
        r"""
        Return a function that returns the object given a FindStat
        representation.

        OUTPUT:

        The function that produces the Sage object given its FindStat
        representation as a string.

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: c = FindStatCollection("Posets")                              # optional -- internet
            sage: p = c.from_string()('([(0, 2), (2, 1)], 3)')                  # optional -- internet
            sage: p.cover_relations()                                           # optional -- internet
            [[0, 2], [2, 1]]

            sage: c = FindStatCollection("Binary Words")                        # optional -- internet
            sage: w = c.from_string()('010101')                                 # optional -- internet
            sage: w in c._data["Code"].elements_on_level(6)                     # optional -- internet
            True
        """
        return self._data["Code"].string_to_element

    def _repr_(self):
        r"""
        Return the representation of the FindStat collection.

        OUTPUT: the representation, including the identifier and the name

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Binary trees")                            # optional -- internet
            Cc0010: Binary trees
        """
        if self._sageconstructor_overridden:
            return "a subset of %s: %s" % (self.id_str(), self._data["NamePlural"])
        return "%s: %s" % (self.id_str(), self._data["NamePlural"])

    def name(self, style='singular'):
        r"""
        Return the name of the FindStat collection.

        INPUT:

        - ``style`` -- string (default: ``'singular'``); can be
          ``'singular'``, or ``'plural'``

        OUTPUT: the name of the FindStat collection, in singular or in plural

        EXAMPLES::

            sage: from sage.databases.findstat import FindStatCollection
            sage: FindStatCollection("Binary trees").name()                     # optional -- internet
            'Binary tree'

            sage: FindStatCollection("Binary trees").name(style='plural')       # optional -- internet
            'Binary trees'
        """
        if style == "singular":
            return self._data["Name"]
        if style == "plural":
            return self._data["NamePlural"]
        raise ValueError("argument 'style' (=%s) must be 'singular' or 'plural'" % style)


from collections import namedtuple
_SupportedFindStatCollection = namedtuple("SupportedFindStatCollection",
                                          ["string_to_element",
                                           "element_to_string",
                                           "elements_on_level", # return all elements on given level
                                           "element_level",     # return level of a given element
                                           "is_element"])       # return whether element is member of this collection

# this dictionary must be sorted so that subclasses come before
# superclasses, eg., "StandardTableaux" before "SemistandardTableaux"
_SupportedFindStatCollections = {
    "Permutations":
    _SupportedFindStatCollection(lambda x: Permutation(literal_eval(x)),
                                 str,
                                 Permutations,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, Permutation)),
    "BinaryWords":
    _SupportedFindStatCollection(lambda x: Word((int(e) for e in str(x)), alphabet=[0,1]),
                                 str,
                                 lambda x: Words([0,1], length=x),
                                 lambda x: x.length(),
                                 lambda x: isinstance(x, Word_class)),
    "AlternatingSignMatrices":
    _SupportedFindStatCollection(lambda x: AlternatingSignMatrix(literal_eval(x)),
                                 lambda x: str(list(map(list, x.to_matrix().rows()))),
                                 AlternatingSignMatrices,
                                 lambda x: x.to_matrix().nrows(),
                                 lambda x: isinstance(x, AlternatingSignMatrix)),
    "BinaryTrees":
    _SupportedFindStatCollection(lambda x: BinaryTree(str(x)),
                                 str,
                                 BinaryTrees,
                                 lambda x: x.node_number(),
                                 lambda x: isinstance(x, BinaryTree)),
    "Cores":
    _SupportedFindStatCollection(lambda x: Core(*literal_eval(x)),
                                 lambda X: "( " + X._repr_() + ", " + str(X.k()) + " )",
                                 lambda x: Cores(x[1], x[0]),
                                 lambda x: (x.length(), x.k()),
                                 lambda x: isinstance(x, Core)),
    "DyckPaths":
    _SupportedFindStatCollection(lambda x: DyckWord(literal_eval(x)),
                                 lambda x: str(list(DyckWord(x))),
                                 DyckWords,
                                 lambda x: x.semilength(),
                                 lambda x: isinstance(x, DyckWord)),
    "FiniteCartanTypes":
    _SupportedFindStatCollection(lambda x: CartanType(*literal_eval(str(x))),
                                 str,
                                 _finite_irreducible_cartan_types_by_rank,
                                 lambda x: x.rank(),
                                 lambda x: isinstance(x, CartanType_abstract)),
    "GelfandTsetlinPatterns":
    _SupportedFindStatCollection(lambda x: GelfandTsetlinPattern(literal_eval(x)),
                                 str,
                                 lambda x: (P
                                            for la in Partitions(x[1], max_length=x[0])
                                            for P in GelfandTsetlinPatterns(top_row=la + [0]*(x[0]-len(la)))),
                                 lambda x: (len(x[0]), sum(x[0])),
                                 lambda x: (x == GelfandTsetlinPatterns
                                            or isinstance(x, GelfandTsetlinPattern))),
    "Graphs":
    _SupportedFindStatCollection(lambda x: (lambda E, V: Graph([list(range(V)),
                                                                lambda i,j: (i,j) in E or (j,i) in E],
                                                               immutable=True))(*literal_eval(x)),
                                 lambda X: str((X.edges(labels=False, sort=True), X.num_verts())),
                                 lambda x: (g.copy(immutable=True) for g in graphs(x, copy=False)),
                                 lambda x: x.num_verts(),
                                 lambda x: isinstance(x, Graph)),
    "IntegerPartitions":
    _SupportedFindStatCollection(lambda x: Partition(literal_eval(x)),
                                 str,
                                 Partitions,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, Partition)),
    "IntegerCompositions":
    _SupportedFindStatCollection(lambda x: Composition(literal_eval(x)),
                                 str,
                                 Compositions,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, Composition)),
    "OrderedTrees":
    _SupportedFindStatCollection(lambda x: OrderedTree(literal_eval(x)),
                                 str,
                                 OrderedTrees,
                                 lambda x: x.node_number(),
                                 lambda x: isinstance(x, OrderedTree)),
    "ParkingFunctions":
    _SupportedFindStatCollection(lambda x: ParkingFunction(literal_eval(x)),
                                 str,
                                 ParkingFunctions,
                                 len,
                                 lambda x: isinstance(x, ParkingFunction)),
    "Lattices":
    _SupportedFindStatCollection(lambda x: (lambda R, E: LatticePoset((list(range(E)), R)))(*literal_eval(x)),
                                 lambda X: str((sorted(X._hasse_diagram.cover_relations()),
                                                len(X._hasse_diagram.vertices(sort=False)))),
                                 _finite_lattices,
                                 lambda x: x.cardinality(),
                                 lambda x: isinstance(x, FiniteLatticePoset)),
    "Posets":
    _SupportedFindStatCollection(lambda x: (lambda R, E: Poset((list(range(E)), R)))(*literal_eval(x)),
                                 lambda X: str((sorted(X._hasse_diagram.cover_relations()),
                                                len(X._hasse_diagram.vertices(sort=False)))),
                                 Posets,
                                 lambda x: x.cardinality(),
                                 lambda x: isinstance(x, FinitePoset)),
    "StandardTableaux":
    _SupportedFindStatCollection(lambda x: StandardTableau(literal_eval(x)),
                                 str,
                                 StandardTableaux,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, StandardTableau)),
    "SemistandardTableaux":
    _SupportedFindStatCollection(lambda x: SemistandardTableau(literal_eval(x)),
                                 str,
                                 lambda x: (T for T in SemistandardTableaux(size=x[0], max_entry=x[1])
                                            if max(T.entries()) == x[1]),
                                 lambda x: (x.size(), max(x.entries())),
                                 lambda x: isinstance(x, SemistandardTableau)),
    "PerfectMatchings":
    _SupportedFindStatCollection(lambda x: PerfectMatching(literal_eval(x)),
                                 str,
                                 PerfectMatchings,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, PerfectMatching)),
    "SetPartitions":
    _SupportedFindStatCollection(lambda x: SetPartition(literal_eval(x.replace('{','[').replace('}',']'))),
                                 str,
                                 SetPartitions,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, SetPartition)),
    "SkewPartitions":
    _SupportedFindStatCollection(lambda x: SkewPartition(literal_eval(x)),
                                 str,
                                 SkewPartitions,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, SkewPartition)),
    "SignedPermutations":
    _SupportedFindStatCollection(lambda x: SignedPermutations(len(literal_eval(x)))(list(literal_eval(x))),
                                 str,
                                 SignedPermutations,
                                 lambda x: len(list(x)),
                                 lambda x: isinstance(x, SignedPermutation)),
    "PlanePartitions":
    _SupportedFindStatCollection(lambda x: PlanePartition(literal_eval(x)),
                                 lambda X: str(list(X)).replace(" ",""),
                                 _plane_partitions_by_size,
                                 lambda x: sum(sum(la) for la in x),
                                 lambda x: isinstance(x, PlanePartition)),
    "DecoratedPermutations":
    _SupportedFindStatCollection(lambda x: DecoratedPermutation([v if v > 0 else (i if v == 0 else -i)
                                                                 for i, v in enumerate(literal_eval(x.replace("+","0").replace("-","-1")), 1)]),
                                 lambda x: "[" + ",".join((str(v) if abs(v) != i else ("+" if v > 0 else "-")
                                                           for i, v in enumerate(x, 1))) + "]",
                                 DecoratedPermutations,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, DecoratedPermutation)),
    "OrderedSetPartitions":
    _SupportedFindStatCollection(lambda x: OrderedSetPartition(literal_eval(x.replace('{','[').replace('}',']'))),
                                 str,
                                 OrderedSetPartitions,
                                 lambda x: x.size(),
                                 lambda x: isinstance(x, OrderedSetPartition))}


class FindStatCollections(UniqueRepresentation, Parent):
    r"""
    The class of FindStat collections.

    The elements of this class are combinatorial collections in
    FindStat as of January 2020.  If a new collection was added to the
    web service since then, the dictionary ``_SupportedFindStatCollections``
    in this module has to be updated accordingly.

    EXAMPLES::

        sage: from sage.databases.findstat import FindStatCollections
        sage: sorted(c for c in FindStatCollections() if c.is_supported())      # optional -- internet
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
         Cc0023: Parking functions,
         Cc0024: Binary words,
         Cc0025: Plane partitions,
         Cc0026: Decorated permutations,
         Cc0027: Signed permutations,
         Cc0028: Skew partitions,
         Cc0029: Lattices,
         Cc0030: Ordered set partitions]
    """
    def __init__(self):
        """
        Fetch the collections from FindStat.

        TESTS::

            sage: from sage.databases.findstat import FindStatCollections
            sage: C = FindStatCollections()                                     # optional -- internet
            sage: TestSuite(C).run()                                            # optional -- internet
        """
        fields = "LevelsWithSizes,Name,NamePlural,NameWiki"
        url = FINDSTAT_API_COLLECTIONS + "?fields=" + fields
        d = _get_json(url, object_pairs_hook=dict)["included"]["Collections"]
        for id, data in d.items():
            data["LevelsWithSizes"] = {literal_eval(level): size
                                       for level, size in data["LevelsWithSizes"].items()}
            if data["NameWiki"] in _SupportedFindStatCollections:
                data["Code"] = _SupportedFindStatCollections[data["NameWiki"]]
            else:
                print("%s provides a new collection:" % FindStat())
                print("    %s: %s" % (id, data["NamePlural"]))
                print("To use it with this interface, it has to be added to the dictionary")
                print("    _SupportedFindStatCollections in src/sage/databases/findstat.py")
                print("of the SageMath distribution.  Please open an issue on github!")
#                print("Very likely, the following code would work:")
#                fields = "SageCodeElementToString,SageCodeElementsOnLevel,SageCodeStringToElement"
#                url = FINDSTAT_API_COLLECTIONS + id + "?fields=" + fields
#                print(json.load(urlopen(url))["included"]["Collections"][id])

        def position(item):
            try:
                return tuple(_SupportedFindStatCollections).index(item[1]["NameWiki"])
            except ValueError:
                return len(_SupportedFindStatCollections)

        self._findstat_collections = dict(sorted(d.items(), key=position))
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, entry):
        """
        Initialize a FindStat collection.

        INPUT:

        See :class:`FindStatCollection`.

        TESTS:

        Create an object and find its collection::

            sage: from sage.databases.findstat import FindStatCollection, FindStatCollections
            sage: sorted([FindStatCollection(c.first_terms(lambda x: 0)[0][0]) for c in FindStatCollections() if c.is_supported()]) # optional -- internet
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
             Cc0023: Parking functions,
             Cc0024: Binary words,
             Cc0025: Plane partitions,
             Cc0026: Decorated permutations,
             Cc0027: Signed permutations,
             Cc0028: Skew partitions,
             Cc0029: Lattices,
             Cc0030: Ordered set partitions]

            sage: FindStatCollection(Permutation([1,2,3]))                      # optional -- internet
            Cc0001: Permutations

            sage: FindStatCollection(Permutations(3))                           # optional -- internet
            a subset of Cc0001: Permutations

            sage: FindStatCollection(PerfectMatching([[1,2]]))                  # optional -- internet
            Cc0012: Perfect matchings

            sage: FindStatCollection(PerfectMatchings(4))                       # optional -- internet
            a subset of Cc0012: Perfect matchings

            sage: FindStatCollection(SetPartition([[1,2]]))                     # optional -- internet
            Cc0009: Set partitions

            sage: FindStatCollection(CartanType("A3"))                          # optional -- internet
            Cc0022: Finite Cartan types

            sage: cc = FindStatCollection(graphs(3)); cc                        # optional -- internet
            a subset of Cc0020: Graphs
            sage: cc.first_terms(lambda x: x.edges(labels=False, sort=True)).list()        # optional -- internet
            [(Graph on 3 vertices, []),
             (Graph on 3 vertices, [(0, 2)]),
             (Graph on 3 vertices, [(0, 2), (1, 2)]),
             (Graph on 3 vertices, [(0, 1), (0, 2), (1, 2)])]

            sage: len(cc.first_terms(lambda x: x.edges(labels=False, sort=False)).list())   # optional -- internet
            4

        Check that we can override the automatic detection::

            sage: l = [(T, len(T)) for n in range(1,5) for T in StandardTableaux(n)]        # optional -- internet
            sage: qu = findstat("Semistandardtableaux", l, depth=1)                         # optional -- internet
            sage: qu[0]                                                                     # optional -- internet
            St000010oMp00077 (quality [100, 100])
        """
        if isinstance(entry, self.Element):
            return entry

        if isinstance(entry, str):
            # find by name in self._findstat_collections (ignoring case and spaces)
            def normalize(e):
                return "".join(e.split()).upper()

            for id, data in self._findstat_collections.items():
                if normalize(entry) in (normalize(id),
                                        normalize(data["NameWiki"]),
                                        normalize(data["NamePlural"]),
                                        normalize(data["Name"])):
                    return self.element_class(self, id, data, None)

        elif isinstance(entry, (int, Integer)):
            # find by id in _findstat_collections
            for id, data in self._findstat_collections.items():
                if entry == int(id[2:]):
                    return self.element_class(self, id, data, None)

        else:
            # find collection given an object or a constructor

            # it would be good to first check whether entry is an
            # element rather than a type - unfortunately, we cannot
            # test with isinstance(_, SageObject), since this is True
            # for CartanType.

            # first check whether the class fits:
            for id, data in self._findstat_collections.items():
                if ("Code" in data
                    and (data["Code"].is_element(entry)
                         # elements_on_level is rarely equal to entry
                         # (it may be a function), but it is
                         # convenient for some types
                         or data["Code"].elements_on_level == entry)):
                    return self.element_class(self, id, data, None)

            # check whether entry is iterable (it's not a string!)
            # producing a subset of objects

            # WARNING: SkewPartition is iterable and produces two
            # elements of Partition

            # WARNING: we have to remember all elements of the
            # generator, because it is used again, for example in
            # self.first_terms
            try:
                entries = lazy_list(entry)
                obj = entries[0]
            except (TypeError, IndexError):
                pass
            else:
                for id, data in self._findstat_collections.items():
                    if "Code" in data and data["Code"].is_element(obj):
                        return self.element_class(self, id, data, entries)

        raise ValueError("could not find FindStat collection for %s" % str(entry))

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
            sage: sorted(FindStatCollections())[0]                              # optional -- internet
            Cc0001: Permutations
        """
        for c in self._findstat_collections:
            yield FindStatCollection(c)

    Element = FindStatCollection
