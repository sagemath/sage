r"""
ISGCI: Information System on Graph Classes and their Inclusions

This module implements an interface to the
`ISGCI <http://www.graphclasses.org/>`_ database in Sage.

This database gathers information on graph classes and their inclusions in each
other. It also contains information on the complexity of several computational
problems.

It is available on the `GraphClasses.org <http://www.graphclasses.org/>`_
website maintained by H.N. de Ridder et al.

How to use it?
--------------

Presently, it is possible to use this database through the variables and methods
present in the :obj:`graph_classes <GraphClasses>` object.
For instance::

    sage: Trees = graph_classes.Tree
    sage: Chordal = graph_classes.Chordal

Inclusions
^^^^^^^^^^

It is then possible to check the inclusion of classes inside of others, if the
information is available in the database::

    sage: Trees <= Chordal
    True

And indeed, trees are chordal graphs.

The ISGCI database is not all-knowing, and so comparing two classes can return
``True``, ``False``, or ``Unknown`` (see the :mod:`documentation of the Unknown
truth value <sage.misc.unknown>`).

An *unknown* answer to ``A <= B`` only means that ISGCI cannot deduce from the
information in its database that ``A`` is a subclass of ``B`` nor that it is
not. For instance, ISGCI does not know at the moment that some chordal graphs
are not trees::

    sage: graph_classes.Chordal <= graph_classes.Tree
    Unknown

Descriptions
^^^^^^^^^^^^

Given a graph class, one can obtain its associated information in the ISGCI
database with the :meth:`~GraphClass.description` method::

    sage: Chordal.description()
    Class of graphs : Chordal
    -------------------------
    id                             :  gc_32
    name                           :  chordal
    ...
    Problems :
    -----------
    3-Colourability                :  Linear
    Clique                         :  Polynomial
    Clique cover                   :  Polynomial
    ...

It is possible to obtain the complete list of the classes stored in ISGCI by
calling the :meth:`~GraphClasses.show_all` method (beware -- long output)::

    sage: graph_classes.show_all()
    id        | name                                     | type                 | smallgraph
    ----------------------------------------------------------------------------------------------------------------------
    gc_309    | $K_4$--minor--free                       | base                 |
    gc_541    | $N^*$                                    | base                 |
    gc_215    | $N^*$--perfect                           | base                 |
    gc_5      | $P_4$--bipartite                         | base                 |
    gc_3      | $P_4$--brittle                           | base                 |
    gc_6      | $P_4$--comparability                     | base                 |
    gc_7      | $P_4$--extendible                        | base                 |
    ...

Until a proper search method is implemented, this lets one find classes which do
not appear in :obj:`graph_classes.* <GraphClasses>`.

To retrieve a class of graph from its ISGCI ID one may use
the :meth:`~GraphClasses.get_class` method::

    sage: GC = graph_classes.get_class("gc_5")
    sage: GC
    $P_4$--bipartite graphs

Recognition of graphs
^^^^^^^^^^^^^^^^^^^^^

The graph classes represented by the ISGCI database can alternatively be used to
access recognition algorithms. For instance, in order to check that a given
graph is a tree one has the following the options ::

    sage: graphs.PathGraph(5) in graph_classes.Tree
    True

or::

    sage: graphs.PathGraph(5).is_tree()
    True

Furthermore, all ISGCI graph classes which are defined by the exclusion of a
finite sequence of induced subgraphs benefit from a generic recognition
algorithm. For instance ::

    sage: g = graphs.PetersenGraph()
    sage: g in graph_classes.ClawFree
    False
    sage: g.line_graph() in graph_classes.ClawFree
    True

Or directly from ISGCI ::

    sage: gc = graph_classes.get_class("gc_441")
    sage: gc
    diamond--free graphs
    sage: graphs.PetersenGraph() in gc
    True

Predefined classes
------------------

:obj:`graph_classes <GraphClasses>` currently predefines the following graph classes

.. list-table::
   :widths: 20 30
   :header-rows: 1

   * - Class
     - Related methods

   * - Apex

     - :meth:`~sage.graphs.graph.Graph.is_apex`,
       :meth:`~sage.graphs.graph.Graph.apex_vertices`

   * - AT_free

     - :meth:`~sage.graphs.graph.Graph.is_asteroidal_triple_free`

   * - Biconnected

     - :meth:`~sage.graphs.graph.Graph.is_biconnected`,
       :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`,
       :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cuts_tree`

   * - BinaryTrees

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.BalancedTree`,
       :meth:`~sage.graphs.graph.Graph.is_tree`

   * - Bipartite

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.BalancedTree`,
       :meth:`~sage.graphs.generic_graph.GenericGraph.is_bipartite`

   * - Block

     - :meth:`~sage.graphs.graph.Graph.is_block_graph`,
       :meth:`~sage.graphs.generic_graph.GenericGraph.blocks_and_cut_vertices`,
       :meth:`~sage.graphs.graph_generators.GraphGenerators.RandomBlockGraph`

   * - Chordal

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_chordal`

   * - Claw-Free
     - :meth:`~sage.graphs.graph_generators.GraphGenerators.ClawGraph`

   * - Comparability
     -

   * - Gallai

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_gallai_tree`

   * - Grid

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.Grid2dGraph`,
       :meth:`~sage.graphs.graph_generators.GraphGenerators.GridGraph`

   * - Interval

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.RandomIntervalGraph`,
       :meth:`~sage.graphs.graph_generators.GraphGenerators.IntervalGraph`,
       :meth:`~sage.graphs.generic_graph.GenericGraph.is_interval`

   * - Line

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.line_graph_forbidden_subgraphs`,
       :meth:`~sage.graphs.graph.Graph.is_line_graph`

   * - Modular

     - :meth:`~sage.graphs.graph.Graph.modular_decomposition`

   * - Outerplanar

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_circular_planar`

   * - Perfect

     - :meth:`~sage.graphs.graph.Graph.is_perfect`

   * - Planar

     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_planar`

   * - Polyhedral

     - :meth:`~sage.graphs.graph.Graph.is_polyhedral`

   * - Split

     - :meth:`~sage.graphs.graph.Graph.is_split`

   * - Tree

     - :meth:`~sage.graphs.graph_generators.GraphGenerators.trees`,
       :meth:`~Graph.is_tree`

   * - UnitDisk
     - :meth:`~sage.graphs.graph_generators.GraphGenerators.IntervalGraph`

   * - UnitInterval
     - :meth:`~sage.graphs.generic_graph.GenericGraph.is_interval`

Sage's view of ISGCI
--------------------

The database is stored by Sage in two ways.

**The classes**: the list of all graph classes and their properties is stored
in a huge dictionary (see :meth:`~sage.graphs.isgci.GraphClasses.classes`).
Below is what Sage knows of ``gc_249``::

    sage: graph_classes.classes()['gc_249']        # random
    {'problem':
        {'Independent set': 'Polynomial',
         'Treewidth': 'Unknown',
         'Weighted independent set': 'Polynomial',
         'Cliquewidth expression': 'NP-complete',
         'Weighted clique': 'Polynomial',
         'Clique cover': 'Unknown',
         'Domination': 'NP-complete',
         'Clique': 'Polynomial',
         'Colourability': 'NP-complete',
         'Cliquewidth': 'Unbounded',
         '3-Colourability': 'NP-complete',
         'Recognition': 'Linear'},
     'type': 'base',
     'id': 'gc_249',
     'name': 'line'}

**The class inclusion digraph**: Sage remembers the class inclusions through
the inclusion digraph (see :meth:`~sage.graphs.isgci.GraphClasses.inclusion_digraph`).
Its nodes are ID of ISGCI classes::

    sage: d = graph_classes.inclusion_digraph()
    sage: d.vertices(sort=True)[-10:]
    ['gc_990', 'gc_991', 'gc_992', 'gc_993', 'gc_994', 'gc_995', 'gc_996', 'gc_997', 'gc_998', 'gc_999']

An arc from ``gc1`` to ``gc2`` means that ``gc1`` is a superclass of ``gc2``.
This being said, not all edges are stored ! To ensure that a given class is
included in another one, we have to check whether there is in the digraph a
``path`` from the first one to the other::

    sage: bip_id = graph_classes.Bipartite._gc_id
    sage: perfect_id = graph_classes.Perfect._gc_id
    sage: d.has_edge(perfect_id, bip_id)
    False
    sage: d.distance(perfect_id, bip_id)
    2

Hence bipartite graphs are perfect graphs. We can see how ISGCI obtains this
result ::

    sage: p = d.shortest_path(perfect_id, bip_id)
    sage: len(p) - 1
    2
    sage: print(p)                  # random
    ['gc_56', 'gc_76', 'gc_69']
    sage: for c in p:
    ....:     print(graph_classes.get_class(c))
    perfect graphs
    ...
    bipartite graphs

What ISGCI knows is that perfect graphs contain unimodular graph which contain
bipartite graphs. Therefore bipartite graphs are perfect !

.. NOTE::

    The inclusion digraph is **NOT ACYCLIC**. Indeed, several entries exist in
    the ISGCI database which represent the same graph class, for instance
    Perfect graphs and Berge graphs::

        sage: graph_classes.inclusion_digraph().is_directed_acyclic()
        False
        sage: Berge = graph_classes.get_class("gc_274"); Berge
        Berge graphs
        sage: Perfect = graph_classes.get_class("gc_56"); Perfect
        perfect graphs
        sage: Berge <= Perfect
        True
        sage: Perfect <= Berge
        True
        sage: Perfect == Berge
        True

Information for developers
--------------------------

* The database is loaded not *so* large, but it is still preferable to
  only load it on demand. This is achieved through the cached methods
  :meth:`~sage.graphs.isgci.GraphClasses.classes` and
  :meth:`~sage.graphs.isgci.GraphClasses.inclusion_digraph`.

* Upon the first access to the database, the information is extracted
  from the XML file and stored in the cache of three methods:

  * ``sage.graphs.isgci._classes`` (dictionary)
  * ``sage.graphs.isgci._inclusions`` (list of dictionaries)
  * ``sage.graphs.isgci._inclusion_digraph`` (DiGraph)

  Note that the digraph is only built if necessary (for instance if
  the user tries to compare two classes).

.. TODO::

    Technical things:

    * Query the database for non-inclusion results so that comparisons can
      return ``False``, and implement strict inclusions.

    * Implement a proper search method for the classes not listed in
      :obj:`graph_classes <GraphClasses>`

      .. SEEALSO:: :func:`sage.graphs.isgci.show_all`.

    * Some of the graph classes appearing in :obj:`graph_classes
      <GraphClasses>` already have a recognition
      algorithm implemented in Sage. It would be so nice to be able to
      write ``g in Trees``, ``g in Perfect``, ``g in Chordal``, ... :-)

    Long-term stuff:

    * Implement simple accessors for all the information in the ISGCI
      database (as can be done from the website)

    * Implement intersection of graph classes

    * Write generic recognition algorithms for specific classes (when a graph
      class is defined by the exclusion of subgraphs, one can write a generic
      algorithm checking the existence of each of the graphs, and this method
      already exists in Sage).

    * Improve the performance of Sage's graph library by letting it take
      advantage of the properties of graph classes. For example,
      :meth:`Graph.independent_set` could use the library to detect that a given
      graph is, say, a tree or a planar graph, and use a specialized algorithm
      for finding an independent set.

AUTHORS:
--------

* H.N. de Ridder et al. (ISGCI database)
* Nathann Cohen (Sage implementation)

Methods
-------
"""

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import CachedRepresentation, UniqueRepresentation
from sage.misc.unknown import Unknown
from sage.features.databases import DatabaseGraphs
from sage.misc.cachefunc import cached_method

import os
import zipfile
from urllib.request import urlopen
from ssl import create_default_context as default_context

# ****************************************************************************
#      Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         https://www.gnu.org/licenses/
# ****************************************************************************

_XML_FILE = "isgci_sage.xml"
_SMALLGRAPHS_FILE = "smallgraphs.txt"


class GraphClass(SageObject, CachedRepresentation):
    r"""
    An instance of this class represents a Graph Class, matching some entry in
    the ISGCI database.

    EXAMPLES:

    Testing the inclusion of two classes::

        sage: Chordal = graph_classes.Chordal
        sage: Trees = graph_classes.Tree
        sage: Trees <= Chordal
        True
        sage: Chordal <= Trees
        Unknown

    TESTS::

        sage: Trees >= Chordal
        Unknown
        sage: Chordal >= Trees
        True
    """
    def __init__(self, name, gc_id, recognition_function=None):
        r"""
        Class constructor.

        INPUT:

        - ``gc_id`` -- the ISGCI class ID

        - ``recognition_function`` -- a function of one argument `g`, which
          return boolean answers to the question : *does ``g`` belong to the
          class represented by ``gc_id`` ?*

        EXAMPLES::

            sage: graph_classes.Chordal  # indirect doctest
            Chordal graphs
        """
        self._name = name
        self._gc_id = gc_id

        if recognition_function is not None:
            self._recognition_function = recognition_function

    def _repr_(self):
        r"""
        Return a short description of the class.

        EXAMPLES::

            sage: graph_classes.Chordal  # indirect doctest
            Chordal graphs
        """
        return self._name + " graphs"

    def __hash__(self):
        r"""
        Return the class' ID hash.

        EXAMPLES::

            sage: hash(graph_classes.Chordal) == hash(graph_classes.Chordal)
            True
        """
        return hash(self._gc_id)

    def __le__(self, other):
        r"""
        <= operator.

        EXAMPLES::

            sage: graph_classes.Chordal <= graph_classes.Tree
            Unknown
        """
        return other >= self

    def __ge__(self, other):
        r"""
        >= operator.

        EXAMPLES::

            sage: graph_classes.Chordal >= graph_classes.Tree
            True
        """
        inclusion_digraph = GraphClasses().inclusion_digraph()
        if inclusion_digraph.shortest_path(self._gc_id, other._gc_id):
            return True
        else:
            return Unknown

    def __eq__(self, other):
        r"""
        == operator.

        EXAMPLES::

            sage: graph_classes.Chordal == graph_classes.Tree
            Unknown
        """
        return self >= other and other >= self

    def __lt__(self, other):
        r"""
        >, !=, and < operators.

        EXAMPLES::

            sage: graph_classes.Chordal > graph_classes.Tree
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: graph_classes.Chordal < graph_classes.Tree
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: graph_classes.Chordal != graph_classes.Tree
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    __gt__ = __ne__ = __lt__

    def forbidden_subgraphs(self):
        r"""
        Return the list of forbidden induced subgraphs defining the class.

        If the graph class is not defined by a *finite* list of forbidden
        induced subgraphs, ``None`` is returned instead.

        EXAMPLES::

            sage: graph_classes.Perfect.forbidden_subgraphs()
            sage: gc = graph_classes.get_class('gc_62')
            sage: gc
            claw--free graphs
            sage: gc.forbidden_subgraphs()
            [Graph on 4 vertices]
            sage: gc.forbidden_subgraphs()[0].is_isomorphic(graphs.ClawGraph())
            True
        """
        classes = GraphClasses().classes()
        gc = classes[self._gc_id]

        if gc.get("type", None) != "forbidden":
            return None

        excluded = gc.get("smallgraph", None)

        if not excluded:
            return None

        if not isinstance(excluded, list):
            excluded = [excluded]

        smallgraphs = GraphClasses().smallgraphs()

        if not all(g in smallgraphs for g in excluded):
            return None

        return [smallgraphs[g] for g in excluded]

    def __contains__(self, g):
        r"""
        Check if ``g`` belongs to the graph class represented by ``self``.

        EXAMPLES::

            sage: graphs.CompleteBipartiteGraph(3,3) in graph_classes.Bipartite
            True
            sage: graphs.CompleteGraph(4) in graph_classes.Chordal
            True
            sage: graphs.CompleteGraph(4) in graph_classes.Comparability
            True
            sage: graphs.CompleteGraph(4) in graph_classes.Interval
            True
            sage: graphs.CompleteGraph(4) in graph_classes.Line
            True
            sage: graphs.CompleteGraph(4) in graph_classes.Perfect
            True
            sage: graphs.CompleteGraph(4) in graph_classes.Planar
            True
            sage: graphs.CompleteGraph(4) in graph_classes.Split
            True
            sage: graphs.PathGraph(4) in graph_classes.Tree
            True
        """
        from sage.graphs.graph import Graph

        if not isinstance(g, Graph):
            return False

        if hasattr(self, "_recognition_function"):
            return self._recognition_function(g)

        excluded = self.forbidden_subgraphs()

        if excluded is None:
            raise NotImplementedError("No recognition algorithm is available "
                                      "for this class.")

        for gg in excluded:
            if g.subgraph_search(gg, induced=True):
                return False

        return True

    def description(self):
        r"""
        Print the information of ISGCI about the current class.

        EXAMPLES::

            sage: graph_classes.Chordal.description()
            Class of graphs : Chordal
            -------------------------
            id                             :  gc_32
            name                           :  chordal
            ...
            Problems :
            -----------
            3-Colourability                :  Linear
            Clique                         :  Polynomial
            Clique cover                   :  Polynomial
            ...
            Recognition                    :  Linear
            ...
        """
        classes = GraphClasses().classes()
        cls = classes[self._gc_id]

        print("Class of graphs : " + self._name)
        print("-" * (len(self._name) + 18))

        for key, value in sorted(cls.items()):
            if value != "" and key != "problem":
                print("{:30} : {}".format(key, value))

        print("\nProblems :")
        print("-" * 11)

        for pbname, data in sorted(cls["problem"].items()):
            if "complexity" in data:
                print("{:30} : {}".format(pbname, data["complexity"]))


class GraphClasses(UniqueRepresentation):
    def get_class(self, id):
        r"""
        Return the class corresponding to the given id in the ISGCI database.

        INPUT:

        - ``id`` -- string; the desired class' ID

        .. SEEALSO::

            :meth:`~sage.graphs.isgci.GraphClasses.show_all`

        EXAMPLES:

        With an existing id::

            sage: Cographs = graph_classes.get_class("gc_151")
            sage: Cographs
            cograph graphs

        With a wrong id::

            sage: graph_classes.get_class(-1)
            Traceback (most recent call last):
            ...
            ValueError: The given class id does not exist in the ISGCI database. Is the db too old ? You can update it with graph_classes.update_db().
        """
        classes = self.classes()
        if id in classes:
            c = classes[id]

            if c.get("name", ""):
                name = c["name"]
            else:
                name = "class " + str(id)

            return GraphClass(name, id)
        else:
            raise ValueError("The given class id does not exist in the ISGCI "
                             "database. Is the db too old ? You can update it "
                             "with graph_classes.update_db().")

    @cached_method
    def classes(self):
        r"""
        Return the graph classes, as a dictionary.

        Upon the first call, this loads the database from the local XML
        file. Subsequent calls are cached.

        EXAMPLES::

            sage: t = graph_classes.classes()
            sage: type(t)
            <... 'dict'>
            sage: sorted(t["gc_151"].keys())
            ['id', 'name',... 'problem',... 'type']
            sage: t["gc_151"]['name']
            'cograph'
            sage: t["gc_151"]['problem']['Clique']
            {'complexity': 'Linear'}
        """
        self._get_ISGCI()
        return self.classes()

    @cached_method
    def inclusions(self):
        r"""
        Return the graph class inclusions.

        OUTPUT: list of dictionaries

        Upon the first call, this loads the database from the local XML file.
        Subsequent calls are cached.

        EXAMPLES::

            sage: t = graph_classes.inclusions()
            sage: type(t)
            <... 'list'>
            sage: t[0]
            {'sub': 'gc_1', 'super': 'gc_2'}
        """
        self._get_ISGCI()
        return self.inclusions()

    @cached_method
    def smallgraphs(self):
        r"""
        Return a dictionary associating a graph to a graph description string.

        Upon the first call, this loads the database from the local XML files.
        Subsequent calls are cached.

        EXAMPLES::

            sage: t = graph_classes.smallgraphs()
            sage: t['2C_4']
            Graph on 8 vertices
            sage: t['2K_3 + e']
            Graph on 6 vertices
            sage: t['fish']
            Graph on 6 vertices
            sage: t['bull']
            Graph on 5 vertices
        """
        self._get_ISGCI()
        return self.smallgraphs()

    @cached_method
    def inclusion_digraph(self):
        r"""
        Return the class inclusion digraph.

        Upon the first call, this loads the database from the local XML file.
        Subsequent calls are cached.

        EXAMPLES::

            sage: g = graph_classes.inclusion_digraph(); g
            Digraph on ... vertices
        """
        classes = self.classes()
        inclusions = self.inclusions()

        from sage.graphs.digraph import DiGraph
        inclusion_digraph = DiGraph()
        inclusion_digraph.add_vertices(classes.keys())

        for edge in inclusions:
            if edge.get("confidence", "") == "unpublished":
                continue
            inclusion_digraph.add_edge(edge['super'], edge['sub'])

        return inclusion_digraph

    def _download_db(self):
        r"""
        Download the current version of the ISGCI db.

        EXAMPLES::

            sage: graph_classes._download_db()  # optional - internet
        """
        import tempfile
        data_dir = os.path.dirname(DatabaseGraphs().absolute_filename())
        u = urlopen('https://www.graphclasses.org/data.zip',
                    context=default_context())
        with tempfile.NamedTemporaryFile(suffix='.zip') as f:
            f.write(u.read())
            z = zipfile.ZipFile(f.name)

            # Save a systemwide updated copy whenever possible
            try:
                z.extract(_XML_FILE, data_dir)
                z.extract(_SMALLGRAPHS_FILE, data_dir)
            except OSError:
                pass

    def _parse_db(self):
        r"""
        Parse the ISGCI database and stores its content in ``self``.

        EXAMPLES::

            sage: graph_classes._parse_db()
        """
        import xml.etree.ElementTree as ET
        from sage.graphs.graph import Graph

        data_dir = os.path.dirname(DatabaseGraphs().absolute_filename())
        xml_file = os.path.join(data_dir, _XML_FILE)
        tree = ET.ElementTree(file=xml_file)
        root = tree.getroot()
        DB = _XML_to_dict(root)

        classes = {c['id']: c for c in DB['GraphClasses']["GraphClass"]}
        for c in classes.values():
            c["problem"] = {pb.pop("name"): pb for pb in c["problem"]}

        inclusions = DB['Inclusions']['incl']

        # Parses the list of ISGCI small graphs
        smallgraph_file = open(os.path.join(data_dir, _SMALLGRAPHS_FILE))
        smallgraphs = {}

        for line in smallgraph_file.readlines():
            key, string = line.split("\t")
            smallgraphs[key] = Graph(string)

        smallgraph_file.close()

        self.inclusions.set_cache(inclusions)
        self.classes.set_cache(classes)
        self.smallgraphs.set_cache(smallgraphs)

    def update_db(self):
        r"""
        Update the ISGCI database by downloading the latest version from
        internet.

        This method downloads the ISGCI database from the website
        `GraphClasses.org <http://www.graphclasses.org/>`_. It then extracts the
        zip file and parses its XML content. The XML file is saved in the directory
        controlled by the :class:`DatabaseGraphs` class (usually, ``$HOME/.sage/db``).

        EXAMPLES::

            sage: graph_classes.update_db()  # optional - internet
            Database downloaded
        """
        self._download_db()

        print("Database downloaded")

        self.classes.clear_cache()
        self.inclusions.clear_cache()
        self.inclusion_digraph.clear_cache()

    def _get_ISGCI(self):
        r"""
        Return the contents of the ISGCI database.

        This method is mostly for internal use, but often provides useful
        information during debugging operations.

        OUTPUT:

        A pair ``(classes, inclusions)`` where ``classes`` is a dict of dict,
        and ``inclusions`` is a list of dicts.

        .. NOTE::

            This method returns the data contained in the most recent ISGCI
            database present on the computer. See :meth:`update_db` to update
            the latter.

        EXAMPLES::

            sage: graph_classes._get_ISGCI()  # long time (4s on sage.math, 2012)
        """
        self._parse_db()

    def show_all(self):
        r"""
        Print all graph classes stored in ISGCI.

        EXAMPLES::

            sage: graph_classes.show_all()
            id        | name                                     | type                 | smallgraph
            ----------------------------------------------------------------------------------------------------------------------
            gc_309    | $K_4$--minor--free                       | base                 |
            gc_541    | $N^*$                                    | base                 |
            gc_215    | $N^*$--perfect                           | base                 |
            gc_5      | $P_4$--bipartite                         | base                 |
            gc_3      | $P_4$--brittle                           | base                 |
            gc_6      | $P_4$--comparability                     | base                 |
            gc_7      | $P_4$--extendible                        | base                 |
            ...
        """
        classes = self.classes()

        # We want to print the different fields, and this dictionary stores the
        # maximal number of characters of each field.
        MAX = {
            "id": 0,
            "type": 0,
            "smallgraph": 0,
            "name": 0
            }

        # We sort the classes alphabetically, though we would like to display
        # the meaningful classes at the top of the list
        def sort_key(x):
            name = x.get("name", "zzzzz")
            return "{}{:4}".format(name, int(x["id"].split('_')[1]))

        classes_list = sorted(classes.values(), key=sort_key)

        # Maximum width of a field
        MAX_LEN = 40

        # Computing the max of each field with the database
        for key in MAX:
            MAX[key] = len(max((str(x.get(key, "")) for x in classes_list), key=len))

        # At most MAX characters per field
        for key, length in MAX.items():
            MAX[key] = min(length, MAX_LEN)

        # Head of the table
        st = ("{:" + str(MAX["id"]) + "}").format("id")
        st += (" | {:" + str(MAX["name"]) + "}").format("name")
        st += (" | {:" + str(MAX["type"]) + "}").format("type")
        st += (" | {:" + str(MAX["smallgraph"]) + "}").format("smallgraph")
        print(st)
        print("-" * (sum(MAX.values())+9))

        # Entries
        for entry in classes_list:
            ID = entry.get("id", "")
            name = entry.get("name", "")
            typ = entry.get("type", "")
            smallgraph = entry.get("smallgraph", "")
            st = ("{:" + str(MAX["id"]) + "}").format(ID)
            st += (" | {:" + str(MAX["name"]) + "}").format(name[:MAX_LEN])
            st += (" | {:" + str(MAX["type"]) + "}").format(typ[:MAX_LEN])
            st += " | " + str(smallgraph)[:MAX_LEN]
            print(st)


def _XML_to_dict(root):
    r"""
    Return the XML data as a dictionary.

    INPUT:

    - ``root`` -- an ``xml.etree.cElementTree.ElementTree`` object

    OUTPUT: a dictionary representing the XML data

    EXAMPLES::

        sage: graph_classes.Perfect.description()  # indirect doctest
        Class of graphs : Perfect
        -------------------------
        id                             :  gc_56
        name                           :  perfect
        ...
    """
    ans = root.attrib.copy()
    for child in root:
        if child.tag in ans:
            if not isinstance(ans[child.tag], list):
                ans[child.tag] = [ans[child.tag]]
            ans[child.tag].append(_XML_to_dict(child))
        else:
            ans[child.tag] = _XML_to_dict(child)

    # If the dictionary is empty, perhaps the only content is a text, and we
    # return this instead. Useful sometimes in the ISGCI db, for graph names.
    if not ans:
        return root.text
    return ans


graph_classes = GraphClasses()

# Any object added to this list should also appear in the class' documentation, at the top of the file.
graph_classes.Apex = GraphClass("Apex", "gc_1181", recognition_function=lambda x: x.is_apex())
graph_classes.AT_free = GraphClass("AT-free", "gc_61", recognition_function=lambda x: x.is_asteroidal_triple_free())
graph_classes.Biconnected = GraphClass("Biconnected", "gc_771", recognition_function=lambda x: x.is_biconnected())
graph_classes.BinaryTrees = GraphClass("BinaryTrees", "gc_847")
graph_classes.Bipartite = GraphClass("Bipartite", "gc_69", recognition_function=lambda x: x.is_bipartite())
graph_classes.Block = GraphClass("Block", "gc_93", recognition_function=lambda x: x.is_block_graph())
graph_classes.Cactus = GraphClass("Cactus", "gc_108", recognition_function=lambda x: x.is_cactus())
graph_classes.Chordal = GraphClass("Chordal", "gc_32", recognition_function=lambda x: x.is_chordal())
graph_classes.ClawFree = GraphClass("Claw-free", "gc_62")
graph_classes.CoGraph = GraphClass("CoGraph", "gc_151", recognition_function=lambda x: x.is_cograph())
graph_classes.Comparability = GraphClass("Comparability", "gc_72", recognition_function=lambda x: x.is_comparability())
graph_classes.DistanceRegular = GraphClass("Distance Regular", "gc_1148", recognition_function=lambda x: x.is_distance_regular())
graph_classes.Gallai = GraphClass("Gallai", "gc_73")
graph_classes.Grid = GraphClass("Grid", "gc_464")
graph_classes.Hamiltonian = GraphClass("Hamiltonian", "gc_1092", recognition_function=lambda x: x.is_hamiltonian())
graph_classes.Interval = GraphClass("Interval", "gc_234", recognition_function=lambda x: x.is_interval())
graph_classes.Line = GraphClass("Line", "gc_249", recognition_function=lambda x: x.is_line_graph())
graph_classes.Modular = GraphClass("Modular", "gc_50")
graph_classes.Outerplanar = GraphClass("Outerplanar", "gc_110")
graph_classes.Perfect = GraphClass("Perfect", "gc_56", recognition_function=lambda x: x.is_perfect())
graph_classes.Permutation = GraphClass("Permutation", "gc_23", recognition_function=lambda x: x.is_permutation())
graph_classes.Planar = GraphClass("Planar", "gc_43", recognition_function=lambda x: x.is_planar())
graph_classes.Polyhedral = GraphClass("Polyhedral", "gc_986", recognition_function=lambda x: x.is_polyhedral())
graph_classes.Split = GraphClass("Split", "gc_39", recognition_function=lambda x: x.is_split())
graph_classes.StronglyRegular = GraphClass("Strongly Regular", "gc_1185", recognition_function=lambda x: x.is_strongly_regular())
graph_classes.Tree = GraphClass("Tree", "gc_342", recognition_function=lambda x: x.is_tree())
graph_classes.TriangleFree = GraphClass("Triangle Free", "gc_371", recognition_function=lambda x: x.is_triangle_free())
graph_classes.UnitDisk = GraphClass("UnitDisk", "gc_389")
graph_classes.UnitInterval = GraphClass("UnitInterval", "gc_299")
graph_classes.WeaklyChordal = GraphClass("Weakly Chordal", "gc_14", recognition_function=lambda x: x.is_weakly_chordal())
