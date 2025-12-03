# -*- coding: utf-8 -*-
r"""
The House of Graphs (HoG)

You can query the HoG (House of Graphs) through Sage in order
to import a graph (identified by its numerical ID in HoG) and
to access its calculated invariants.

AUTHORS:

- Darij Grinberg and Jan Goedgebeur (2023): initial version.

EXAMPLES::

    sage: hog
    The House of Graphs (https://houseofgraphs.org/)

We load the Petersen graph from HoG::

    sage: # optional -- internet
    sage: H = hog(660); H
    660: Petersen Graph
    sage: P = H.to_graph(); P
    Graph on 10 vertices
    sage: P.is_isomorphic(graphs.PetersenGraph())
    True
    sage: H.invariants()[25]
    4
    sage: hog.invariants()[25]
    {'definition': 'The number of <u>edges</u> of the <i>longest induced path</i> of the graph.',
     'invariantId': 25,
     'invariantName': 'Longest Induced Path',
     'keyword': 'LongestInducedPath',
     'typeName': 'i'}
    sage: H.data()["adjacencyList"]
    [[1, 2, 3], [0, 4, 5], [0, 6, 9], [0, 7, 8],
     [1, 6, 8], [1, 7, 9], [2, 4, 7], [3, 5, 6],
     [3, 4, 9], [2, 5, 8]]

A synonym for ``hog(n)`` is ``hog.find_by_id(n)``::

    sage: # optional -- internet
    sage: hog.find_by_id(660)
    660: Petersen Graph

The ID number can also be inputted as a string::

    sage: # optional -- internet
    sage: hog.find_by_id("660")
    660: Petersen Graph
    sage: hog("660")
    660: Petersen Graph

Not all graphs have proper names::

    sage: # optional -- internet
    sage: hog(330)
    330: Graph on 5 vertices

Not all IDs belong to graphs::

    sage: # optional -- internet
    sage: hog(23)
    ... 404 Client Error: ...

The HoG contains a list of known invariants for its graphs.
Some bare-bones support for these invariants is provided in this
module::

    sage: # optional -- internet
    sage: hog.invariants()[5]['invariantName']
    'Clique Number'
    sage: hog.invariants()[5]['definition']
    'The <i>clique number</i> of a graph <i>G</i> is the maximum cardinality of a clique in <i>G</i>.'
    sage: hog(26).invariants()[5]
    4
    sage: hog.invariants()[41]['keyword']
    'Hypohamiltonian'
    sage: hog(26).invariants()[41]
    False

.. TODO: Find by invariants or name etc.

.. TODO: Any other metadata? Comments? Owners?

.. TODO: Should HoGGraph extend Graph? I'm not sure what the
   advantages are; to me it appears that it might confuse the
   reader about what methods come from what class.

Classes and methods
-------------------
"""

# ****************************************************************************
#       Copyright (C) 2023 Darij Grinberg <darijgrinberg!gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import requests

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer import Integer
from sage.matrix.constructor import matrix
from sage.graphs.graph import Graph
from sage.rings.integer_ring import ZZ

hog_url = 'https://houseofgraphs.org/'
hog_api_url = 'https://houseofgraphs.org/api/'


class HoG:
    r"""
    The House of Graphs.

    ``HoG`` is a class representing the House of Graphs. You can
    query it using its methods, but ``HoG`` can also be called
    directly with three arguments:

    - ``query`` -- it can be:

      - a string representing an HoG ID (e.g. '26').
      - an integer representing an HoG ID (e.g. 26).

    OUTPUT:

    - if ``query`` is an integer or a string, returns
      the associated HoG graph.

    EXAMPLES::

        sage: hog
        The House of Graphs (https://houseofgraphs.org/)

    A particular graph can be called by its ID::

        sage: # optional -- internet
        sage: hog('100')
        100: Graph on 6 vertices
        sage: hog(100)
        100: Graph on 6 vertices

    """

    def __call__(self, query, max_results=3, first_result=0):
        r"""
        See the documentation of :class:`HoG`.

        TESTS::

            sage: hog()
            Traceback (most recent call last):
            ...
            TypeError: ...__call__() ...
        """
        return self.find_by_id(query)
        # This will become significantly more complicated
        # if lookup by name or similar is implemented.

    def __repr__(self) -> str:
        r"""
        Return the representation of ``self``.

        TESTS::

            sage: hog
            The House of Graphs (https://houseofgraphs.org/)
        """
        return "The House of Graphs (%s)" % hog_url

    def find_by_id(self, ident, fetch=False):
        r"""

        INPUT:

        - ``ident`` -- the number of the graph in the database (either as
          integer or as string)

        - ``fetch`` -- (bool, default: ``False``) whether to force fetching the
          content of the sequence on the internet.

        OUTPUT:

        - The HoG sequence whose A-number or number corresponds to ``ident``.

        EXAMPLES::

            sage: # optional -- internet
            sage: hog.find_by_id('100')
            100: Graph on 6 vertices
            sage: hog.find_by_id(100)
            100: Graph on 6 vertices
        """
        graph = HoGGraph(ident)
        if fetch:
            graph.online_update()
        return graph

    def invariants(self):
        r"""
        Output a dictionary containing all graph invariants known to
        the HoG database.
        The keys of this dictionary are positive integers.
        The values are graph invariants.
        Each value is itself a dictionary, consisting of its key number
        (key ``invariantId``), a verbose definition (key ``definition``),
        a quick one-word name (key ``keyword``), a longer (mathematically
        meaningful) name (key ``invariantName``), and a type (key
        ``typeName``). The type is one of ``b`` (boolean), ``i``
        (integer) and ``r`` (float).

        EXAMPLES::

            sage: # optional -- internet
            sage: hog.invariants()[43]
            {'definition': 'The number of spanning trees of a graph.
             A <i>spanning tree</i> of a graph is a subgraph that is
             a tree and which contains all vertices of that graph
             and some (or all) of the edges.',
             'invariantId': 43,
             'invariantName': 'Number of Spanning Trees',
             'keyword': 'NumberOfSpanningTrees',
             'typeName': 'i'}
        """
        try:
            return self._invar
        except AttributeError:
            self.online_update(invariants=True)
            return self._invar

    def online_update(self, invariants=False):
        r"""
        Fetch the online HoG to update general data (e.g., a list of
        all invariants).

        TESTS::

            sage: # optional -- internet
            sage: h = hog; h
            The House of Graphs (https://houseofgraphs.org/)
            sage: h.online_update(invariants=True)
            sage: h._invar[13]
            {'definition': 'The <i>domination number</i> of a graph
             <i>G</i> is the minimum cardinality of a set of
             vertices such that every vertices is either in the set
             or is a neighbor of a vertex in the set.',
             'invariantId': 13,
             'invariantName': 'Domination Number',
             'keyword': 'DominationNumber',
             'typeName': 'i'}
        """
        if invariants:
            url = hog_api_url + "invariants"
            js_response = requests.get(url)
            raw = js_response.json()
            rawlist = raw["_embedded"]["invariantModelList"]
            d = {int(rawi["entity"]["invariantId"]): rawi["entity"]
                 for rawi in rawlist}
            self._invar = d

    # The types an invariant can have:
    invariant_types = {"b": bool, "i": int, "r": float}

    def browse(self):
        r"""
        Open the HoG web page in a browser.

        EXAMPLES::

            sage: hog.browse()                     # optional -- internet webbrowser
        """
        import webbrowser
        webbrowser.open(hog_url)


class HoGGraph(SageObject, UniqueRepresentation):
    r"""
    The class of HoG graphs.

    This class implements HoG graphs. They are usually produced by calls to
    the House of Graphs, represented by the class :class:`HoG`.

    .. NOTE::

        A HoG graph (i.e., an instance of this class) is
        not the same as a SageMath graph. To turn the
        former into the latter, use the :meth:`to_graph`
        method on the former.

    .. automethod:: __call__
    """
    @staticmethod
    def __classcall__(cls, ident):
        r"""
        Canonicalize the ID of the graph into a string.

        TESTS::

            sage: hog(45) is hog('45')
            True
        """
        if not isinstance(ident, str):
            ident = str(ident)
        return super().__classcall__(cls, ident)

    def __init__(self, ident):
        r"""
        Initialize a HoG graph.

        There is no fetching of additional information about the graph at
        this point; only the ID number is required to construct a graph.
        In particular, the number might be invalid (i.e., not correspond
        to any item in the HoG), which will raise exceptions once any
        nontrivial methods are called on the result.

        INPUT:

        - ``ident`` -- a string representing the number of the graph or an
          integer representing its number.

        TESTS::

            sage: from sage.databases.hog import HoGGraph
            sage: HoGGraph(26)
            26: Graph on 7 vertices
            sage: HoGGraph("26")
            26: Graph on 7 vertices
            sage: HoGGraph(26) is HoGGraph("26")
            True
        """
        self._id = ident

    def online_update(self, target="data"):
        r"""
        Fetch the online HoG to update the informations about this graph.

        TESTS::

            sage: s = hog("26")
            sage: if hasattr(s, "_data"): delattr(s, "_data")
            sage: s._data
            Traceback (most recent call last):
            ...
            AttributeError: 'HoGGraph' object has no attribute '_data'
            sage: s.online_update(target="data")
            sage: s._data
            {...
             'adjacencyList': [[1, 2, 5, 6],
              [0, 3, 5, 6],
              [0, 4, 5, 6],
              [1, 4, 5, 6],
              [2, 3, 5, 6],
              [0, 1, 2, 3, 4, 6],
              [0, 1, 2, 3, 4, 5]],
             ...}

            sage: if hasattr(s, "_invar"): delattr(s, "_invar")
            sage: s._invar
            Traceback (most recent call last):
            ...
            AttributeError: 'HoGGraph' object has no attribute '_invar'
            sage: s.online_update(target="invariants")
            sage: s._invar[27]
            15

            sage: if hasattr(s, "_embeds"): delattr(s, "_embeds")
            sage: s.online_update(target="embeddings")
            sage: s._embeds[1]
            [...]
        """
        if target == "data":
            url = hog_api_url + "graphs/" + self._id
            js_response = requests.get(url)
            js_response.raise_for_status()
            self._data = js_response.json()

        elif target == "invariants":
            # Fetch the invariants of self.
            url = hog_api_url + "graphs/" + self._id + "/invariants"
            js_response = requests.get(url)
            js_response.raise_for_status()
            ll = js_response.json()['_embedded']['graphInvariantModelList']

            # Cast the invariants into the appropriate types.
            def invar_from_entry(i):
                ident = int(i["entity"]["invariantId"])
                itype = hog.invariant_types[hog.invariants()[ident]["typeName"]]
                return itype(i["entity"]["invariantValue"])
            self._invar = {int(i["entity"]["invariantId"]):
                           invar_from_entry(i)
                           for i in ll if i["entity"]["invariantStatus"] == 2}
            # invariantStatus == 2 means that the invariant has
            # been computed on the HoS server. Other statuses
            # mean uncomputed invariants.

        elif target == "embeddings":
            url = hog_api_url + "graphs/" + self._id + "/embedding"
            js_response = requests.get(url)
            js_response.raise_for_status()
            ll = js_response.json()['_embedded']['embeddingModelList']
            edict = {emb["entity"]["embeddingId"]: emb["entity"]["embedding"]
                     for emb in ll}
            self._embeds = edict

    def to_graph(self, embedding=None):
        r"""
        Return the actual graph (a :cls:``Graph`` object) underlying
        the HoG graph ``self``.

        This forgets all the information downloaded from HoG except
        the graph itself and its standard embedding (in the
        2-dimensional plane). Optionally, a different embedding can
        be provided through the optional ``embedding`` parameter.

        EXAMPLES::

            sage: # optional -- internet
            sage: s = hog("660")
            sage: s.to_graph()
            Graph on 10 vertices
            sage: s.all_embeddings()
            {...}
            sage: s.to_graph(embedding=21566)
            Graph on 10 vertices
        """
        AM = self.data()["adjacencyMatrix"]
        AM = matrix(ZZ, [[int(i) for i in r] for r in AM])
        if embedding is False:
            return Graph(AM, format="adjacency_matrix")
        return Graph(AM, format="adjacency_matrix", pos=self.embedding(embedding=embedding))

    def invariants(self):
        r"""
        Return the graph invariants of HoG graph ``self``,
        retrieved from the HoG database.

        Use ``hog.invariants()`` for the meanings of these invariants.

        EXAMPLES:

        Let us compute the domination number (HoG invariant 13)
        of the Petersen graph (HoG graph 660)::

            sage: # optional -- internet
            sage: hog.invariants()[13]['invariantName']
            'Domination Number'
            sage: hog(660).invariants()[13]
            3
        """
        try:
            return self._invar
        except AttributeError:
            self.online_update(target="invariants")
            return self._invar

    def embedding(self, embedding=None):
        r"""
        Return the standard embedding of the HoG graph ``self``
        given in the HoG database. This is a dictionary of points in
        the plane (provided as pairs of coordinates), indexed by
        vertices of ``self``.

        Note that HoG might have several embeddings of ``self`` in
        the database; they can all be retrieved using the
        :meth:`all_embeddings` method.

        Of course, these embeddings are not guaranteed to be
        planar, even when the graph is planar.

        EXAMPLES::

            sage: # optional -- internet
            sage: hog(660).embedding()
            {0: [-0.5053758991724633, -0.46544052312593087],
             1: [0.13454805435491934, 0.09885219713636184],
             2: [-0.9012033548216514, -0.528004241038393],
             3: [0.30775281428610235, -0.4610066018006447],
             4: [0.39093485483467205, 0.4074994557021103],
             5: [-0.09882598543538978, -0.7277049782582388],
             6: [-0.5736315646159965, 0.3983276914207079],
             7: [-0.3358858825684776, 0.08266786853513786],
             8: [0.7429867105026631, -0.5193210275035198],
             9: [-0.09826837570345093, -1.0987353451584003]}
        """
        if embedding is None:
            M = self.data()["embedding"]
        else:
            M = self.all_embeddings()[embedding]
        return {i: Mi for i, Mi in enumerate(M)}

    def all_embeddings(self):
        r"""
        Return a dictionary containing all embeddings of ``self``
        known to the HoG database.

        This includes the standard embedding (:meth:``embedding``).

        EXAMPLES::

            sage: # optional -- internet
            sage: hog(660).all_embeddings()[21566]
            [[-0.8775468891588245, 0.2732669045189495],
             [-0.2042942155151185, 0.7677182601798132],
             [-0.8945791887222087, -0.511324531389524],
             [0.11697358590230778, 0.26341791868239883],
             [-0.6804465206488972, -0.12330468950606965],
             [0.6935614198570197, 0.2629106524019438],
             [-0.15835548917128717, -1.042450003238761],
             [0.6864200745668088, -0.536005156425335],
             [-0.1961968555594359, -0.12129524244822742],
             [0.12418193748813544, -0.5295397427681126]]
        """
        try:
            return self._embeds
        except AttributeError:
            self.online_update(target="embeddings")
            return self._embeds

    def id(self, format='str'):
        r"""
        The ID of the sequence ``self`` is the A-number that identifies
        ``self``.

        INPUT:

        - ``format`` -- (string, default: 'str').

        OUTPUT:

        - if ``format`` is set to 'str', returns a string.
        - if ``format`` is set to 'int' returns an integer.

        EXAMPLES::

            sage: # optional -- internet
            sage: f = hog(660)
            sage: f.id()
            '660'
            sage: f.id("int")
            660
        """
        if format == 'str':
            return self._id
        elif format == 'int':
            return int(self._id)

    def __hash__(self):
        r"""
        Return the hash of ``self``, which is its numerical HoG ID.

        This method allows unique representation of HoG graphs.

        OUTPUT:

        - Python `int`.

        EXAMPLES::

            sage: # optional -- internet
            sage: s = hog(660)
            sage: hash(s)
            660

        We have unique representation::

            sage: # optional -- internet
            sage: t = hog(660)
            sage: s is t
            True
            sage: s == t
            True
        """
        return self.id(format='int')

    def data(self):
        r"""
        Return the data of the graph ``self``, as a dictionary
        structured in the HoG format.

        The raw entry is fetched online if needed.

        OUTPUT:

        - dictionary.

        EXAMPLES::

            sage: # optional -- internet
            sage: s = hog(660); s
            660: Petersen Graph
            sage: Pet = hog(660)
            sage: Pet.data()
            {...}
            sage: Pet.data()["entity"]["graphId"]
            660
            sage: Pet.data()["entity"]["graphName"]
            'Petersen Graph'
            sage: Pet.data()["entity"]["canonicalForm"]
            'SXNQQE9rV0hH'
            sage: Pet.data()["adjacencyMatrix"]
            [[False, True, True, True, False, False, False, False, False, False],
              [True, False, False, False, True, True, False, False, False, False],
              [True, False, False, False, False, False, True, False, False, True],
              [True, False, False, False, False, False, False, True, True, False],
              [False, True, False, False, False, False, True, False, True, False],
              [False, True, False, False, False, False, False, True, False, True],
              [False, False, True, False, True, False, False, True, False, False],
              [False, False, False, True, False, True, True, False, False, False],
              [False, False, False, True, True, False, False, False, False, True],
              [False, False, True, False, False, True, False, False, True, False]]
            sage: Pet.data()["adjacencyList"]
            [[1, 2, 3],
             [0, 4, 5],
             [0, 6, 9],
             [0, 7, 8],
             [1, 6, 8],
             [1, 7, 9],
             [2, 4, 7],
             [3, 5, 6],
             [3, 4, 9],
             [2, 5, 8]]
        """
        try:
            return self._data
        except AttributeError:
            self.online_update(target="data")
            return self._data

    def name(self) -> str:
        r"""
        Return the name of the graph ``self``.

        This is not a unique identifier.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: # optional -- internet
            sage: hog(660).name()
            'Petersen Graph'
            sage: hog(26).name()
            'Graph on 7 vertices'
        """
        sde = self.data()["entity"]
        if "graphName" in sde:
            if sde["graphName"] is not None:
                return sde["graphName"]
        return repr(self.to_graph())

    def _repr_(self):
        r"""
        Print the name of the graph and its number in HoG.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: # optional -- internet
            sage: hog(660)
            660: Petersen Graph
            sage: hog(26)
            26: Graph on 7 vertices
        """
        return "%s: %s" % (self.id(), self.name())

    def url(self):
        r"""
        Return the URL of the page associated to the graph ``self``.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: # optional -- internet
            sage: hog(26).url()
            'https://houseofgraphs.org/graphs/26'
            sage: hog(55).url()
            'https://houseofgraphs.org/graphs/55'
        """
        return hog_url + "graphs/" + self.id()

    def browse(self):
        r"""
        Open the HoG web page associated to the graph ``self`` in a browser.

        EXAMPLES::

            sage: hog(45).browse()                     # optional -- internet webbrowser

        """
        import webbrowser
        webbrowser.open(self.url())


hog = HoG()
