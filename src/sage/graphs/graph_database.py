"""
Graph database

This module implements classes (:class:`~GraphDatabase`, :class:`~GraphQuery`,
:class:`~GenericGraphQuery`) for interfacing with the ``sqlite`` database
``graphs.db``.

The :class:`~GraphDatabase` class interfaces with the ``sqlite`` database
``graphs.db``. It is an immutable database that inherits from
:class:`~sage.databases.sql_db.SQLDatabase` (see :mod:`sage.databases.sql_db`).

The database contains all unlabeled graphs with 7 or fewer nodes.  This class
will also interface with the optional database package containing all unlabeled
graphs with 8 or fewer nodes. The database(s) consists of five tables, and has
the structure given by the function :func:`~graph_db_info` (For a full
description including column data types, create a :class:`~GraphDatabase`
instance and call the method
:meth:`~sage.databases.sql_db.SQLDatabase.get_skeleton`).

AUTHORS:

- Emily A. Kirkman (2008-09-20): first version of interactive queries,
  cleaned up code and generalized many elements to
  sage.databases.sql_db.py

- Emily A. Kirkman (2007-07-23): inherits GenericSQLDatabase, also
  added classes: GraphQuery and GenericGraphQuery

- Emily A. Kirkman (2007-05-11): initial sqlite version

- Emily A. Kirkman (2007-02-13): initial version (non-sqlite)

REFERENCES:

- Data provided by Jason Grout (Brigham Young University). [Online]
  Available: http://artsci.drake.edu/grout/graphs/
"""

# ##############################################################################
#           Copyright (C) 2007 Emily A. Kirkman
#
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         https://www.gnu.org/licenses/
# ##############################################################################

from . import graph
import os
import re
from sage.rings.integer import Integer
from sage.databases.sql_db import SQLDatabase, SQLQuery
from sage.features.databases import DatabaseGraphs
from sage.graphs.graph import Graph
dblocation = DatabaseGraphs().absolute_filename()


def degseq_to_data(degree_sequence):
    """
    Convert a degree sequence list to a sorted (max-min) integer data type.

    The input degree sequence list (of Integers) is converted to a sorted
    (max-min) integer data type, as used for faster access in the underlying
    database.

    INPUT:

    - ``degree_sequence`` -- list of integers; input degree sequence list

    EXAMPLES::

        sage: from sage.graphs.graph_database import degseq_to_data
        sage: degseq_to_data([2,2,3,1])
        3221
    """
    degree_sequence.sort()
    return sum(di * 10**i for i, di in enumerate(degree_sequence))


def data_to_degseq(data, graph6=None):
    """
    Convert a database integer data type to a degree sequence list.

    INPUT:

    - ``data`` -- integer data type (one digit per vertex representing its
      degree, sorted high to low) to be converted to a degree sequence list

    - ``graph6`` -- string (default: ``None``); the ``graph6`` identifier is
      required for all graphs with no edges, so that the correct number of zeros
      is returned.

    EXAMPLES::

        sage: from sage.graphs.graph_database import data_to_degseq
        sage: data_to_degseq(3221)
        [1, 2, 2, 3]
        sage: data_to_degseq(0, 'D??')
        [0, 0, 0, 0, 0]
    """
    degseq = Integer(data).digits(10)
    if not degseq:
        # compute number of 0s in list from graph6 string
        from sage.graphs.generic_graph_pyx import length_and_string_from_graph6
        return length_and_string_from_graph6(str(graph6))[0] * [0]
    return degseq


def graph6_to_plot(graph6):
    """
    Return a ``Graphics`` object from a ``graph6`` string.

    This method constructs a graph from a ``graph6`` string and returns a
    :class:`sage.plot.graphics.Graphics` object with arguments preset for the
    :meth:`sage.plot.graphics.Graphics.show` method.

    INPUT:

    - ``graph6`` -- a ``graph6`` string

    EXAMPLES::

        sage: from sage.graphs.graph_database import graph6_to_plot
        sage: type(graph6_to_plot('D??'))                                               # needs sage.plot
        <class 'sage.plot.graphics.Graphics'>
    """
    g = Graph(str(graph6))
    return g.plot(layout='circular', vertex_size=30, vertex_labels=False, graph_border=False)


def subgraphs_to_query(subgraphs, db):
    """
    Return a GraphQuery object required for the induced_subgraphs parameter.

    This method constructs and returns a :class:`~GraphQuery` object respecting
    the special input required for the ``induced_subgraphs`` parameter.

    INPUT:

    - ``subgraphs`` -- list of strings; the list should be of one of the
      following two formats:

      - ``['one_of', String, ..., String]`` -- will search for graphs containing
        a subgraph isomorphic to *any* of the ``graph6`` strings in the list

      - ``['all_of', String, ..., String]`` -- will search for graphs containing
        a subgraph isomorphic to *each* of the ``graph6`` strings in the list

    - ``db`` -- a :class:`~GraphDatabase`

    .. NOTE::

        This is a helper method called by the :class:`~GraphQuery` constructor
        to handle this special format. This method should not be used on its own
        because it doesn't set any display columns in the query string, causing
        a failure to fetch the data when run.

    EXAMPLES::

        sage: from sage.graphs.graph_database import subgraphs_to_query
        sage: gd = GraphDatabase()
        sage: q = subgraphs_to_query(['all_of', 'A?', 'B?', 'C?'], gd)
        sage: q.get_query_string()
        'SELECT ,,,,,  FROM misc WHERE ( ( misc.induced_subgraphs regexp ? ) AND (
        misc.induced_subgraphs regexp ? ) ) AND ( misc.induced_subgraphs regexp ? )'
    """
    q = GraphQuery(graph_db=db, induced_subgraphs=subgraphs[1])
    if subgraphs[0] == 'all_of':
        for i in range(2, len(subgraphs)):
            q.intersect(GraphQuery(graph_db=db, induced_subgraphs=subgraphs[i]),
                        in_place=True)
    elif subgraphs[0] == 'one_of':
        for i in range(2, len(subgraphs)):
            q.union(GraphQuery(graph_db=db, induced_subgraphs=subgraphs[i]),
                    in_place=True)
    else:
        raise KeyError('unable to initiate query: illegal input format for induced_subgraphs')
    return q


# tables     columns                    input data type     sqlite data type
# -----------------------------------------------------------------------------
aut_grp = ['aut_grp_size',             # Integer           INTEGER
           'num_orbits',               # Integer           INTEGER
           'num_fixed_points',         # Integer           INTEGER
           'vertex_transitive',        # bool              BOOLEAN
           'edge_transitive']          # bool              BOOLEAN
degrees = ['degree_sequence',          # list              INTEGER (see degseq_to_data module function)
           'min_degree',               # Integer           INTEGER
           'max_degree',               # Integer           INTEGER
           'average_degree',           # Real              REAL
           'degrees_sd',               # Real              REAL
           'regular']                  # bool              BOOLEAN
misc = ['vertex_connectivity',      # Integer           INTEGER
        'edge_connectivity',        # Integer           INTEGER
        'num_components',           # Integer           INTEGER
        'girth',                    # Integer           INTEGER
        'radius',                   # Integer           INTEGER
        'diameter',                 # Integer           INTEGER
        'clique_number',            # Integer           INTEGER
        'independence_number',      # Integer           INTEGER
        'num_cut_vertices',         # Integer           INTEGER
        'min_vertex_cover_size',    # Integer           INTEGER
        'num_spanning_trees',       # Integer           INTEGER
        'induced_subgraphs']        # String            STRING
spectrum = ['spectrum',                 # String            STRING
            'min_eigenvalue',           # Real              REAL
            'max_eigenvalue',           # Real              REAL
            'eigenvalues_sd',           # Real              REAL
            'energy']                   # Real              REAL
graph_data = ['complement_graph6',        # String            STRING
              'eulerian',                 # bool              BOOLEAN
              'graph6',                   # String            STRING
              'lovasz_number',            # Real              REAL
              'num_cycles',               # Integer           INTEGER
              'num_edges',                # Integer           INTEGER
              'num_hamiltonian_cycles',   # Integer           INTEGER
              'num_vertices',             # Integer           INTEGER
              'perfect',                  # bool              BOOLEAN
              'planar']                   # bool              BOOLEAN

valid_kwds = aut_grp + degrees + misc + spectrum + graph_data


def graph_db_info(tablename=None):
    """
    Return a dictionary of allowed table and column names.

    INPUT:

    - ``tablename`` -- restricts the output to a single table

    EXAMPLES::

        sage: sorted(graph_db_info())
        ['aut_grp', 'degrees', 'graph_data', 'misc', 'spectrum']

    ::

        sage: graph_db_info(tablename='graph_data')
        ['complement_graph6',
         'eulerian',
         'graph6',
         'lovasz_number',
         'num_cycles',
         'num_edges',
         'num_hamiltonian_cycles',
         'num_vertices',
         'perfect',
         'planar']
    """
    info = {'graph_data': graph_data,
            'aut_grp': aut_grp,
            'degrees': degrees,
            'misc': misc,
            'spectrum': spectrum}
    if tablename is not None:
        info = info[tablename]
    return info


class GenericGraphQuery(SQLQuery):

    def __init__(self, query_string, database=None, param_tuple=None):
        """
        A query for a :class:`~GraphDatabase`.

        INPUT:

        - ``query_string`` -- string representing the SQL query

        - ``database`` -- (default: ``None``) the :class:`~GraphDatabase`
          instance to query (if ``None`` then a new instance is created)

        - ``param_tuple`` -- tuple of strings (default: ``None``); what to
          replace question marks in ``query_string`` with (optional, but a good
          idea)

        .. NOTE::

           This query class is generally intended for developers and more
           advanced users. It allows you to execute any query, and so may be
           considered unsafe.

        EXAMPLES:

        See :class:`~GraphDatabase` class docstrings or enter::

            sage: G = GraphDatabase()
            sage: G.get_skeleton()
            {...

        to see the underlying structure of the database. Also see
        :class:`sage.databases.sql_db.SQLQuery` in :mod:`sage.databases.sql_db`
        for more info and a tutorial.

        A piece of advice about '?' and param_tuple: it is generally considered
        safer to query with a '?' in place of each value parameter, and using a
        second argument (a tuple of strings) in a call to the ``sqlite``
        database. Successful use of the ``param_tuple`` argument is
        exemplified::

            sage: G = GraphDatabase()
            sage: q = 'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
            sage: param = (22,5)
            sage: Q = SQLQuery(G, q, param)
            sage: Q.show()
            graph_id             graph6               num_vertices         num_edges
            --------------------------------------------------------------------------------
            18                   D??                  5                    0
            19                   D?C                  5                    1
            20                   D?K                  5                    2
            21                   D@O                  5                    2
            22                   D?[                  5                    3
        """
        if database is None:
            database = GraphDatabase()
        if not isinstance(database, GraphDatabase):
            raise TypeError('%s is not a valid GraphDatabase' % database)
        SQLQuery.__init__(self, database, query_string, param_tuple)


class GraphQuery(GenericGraphQuery):

    def __init__(self, graph_db=None, query_dict=None, display_cols=None, **kwds):
        """
        A query for an instance of :class:`~GraphDatabase`.

        This class nicely wraps the :class:`sage.databases.sql_db.SQLQuery`
        class located in :mod:`sage.databases.sql_db` to make the query
        constraints intuitive and with as many pre-definitions as
        possible. (i.e.: since it has to be a :class:`~GraphDatabase`, we
        already know the table structure and types; and since it is immutable,
        we can treat these as a guarantee).

        .. NOTE::

           :class:`sage.databases.sql_db.SQLQuery` functions are available for
           :class:`~GraphQuery`. See :mod:`sage.databases.sql_db` for more
           details.

        INPUT:

        - ``graph_db`` -- :class:`~GraphDatabase` (default: ``None``); instance
          to apply the query to (If ``None``, then a new instance is created)

        - ``query_dict`` -- dictionary (default: ``None``); a dictionary
          specifying the query itself. Format is: ``{'table_name': 'tblname',
          'display_cols': ['col1', 'col2'], 'expression': [col, operator,
          value]}``. If not ``None``, ``query_dict`` will take precedence over
          all other arguments.

        - ``display_cols`` -- list of strings (default: ``None``); a list of
          column names (strings) to display in the result when running or
          showing a query

        - ``kwds`` -- the columns of the database are all keywords. For a
          database table/column structure dictionary, call
          :func:`~graph_db_info`. Keywords accept both single values and lists
          of length 2. The list allows the user to specify an expression other
          than equality. Valid expressions are strings, and for numeric values
          (i.e. Reals and Integers) are: '=','','','=','='. String values also
          accept 'regexp' as an expression argument. The only keyword exception
          to this format is ``induced_subgraphs``, which accepts one of the
          following options:

          - ``['one_of', String, ..., String]`` -- will search for graphs
            containing a subgraph isomorphic to *any* of the ``graph6`` strings
            in the list

          - ``['all_of', String, ..., String]`` -- will search for graphs
            containing a subgraph isomorphic to *each* of the ``graph6`` strings
            in the list

        EXAMPLES::

            sage: Q = GraphQuery(display_cols=['graph6', 'num_vertices', 'degree_sequence'], num_edges=['<=', 5], min_degree=1)
            sage: Q.number_of()
            35
            sage: Q.show()
            Graph6               Num Vertices         Degree Sequence
            ------------------------------------------------------------
            A_                   2                    [1, 1]
            BW                   3                    [1, 1, 2]
            CF                   4                    [1, 1, 1, 3]
            CK                   4                    [1, 1, 1, 1]
            CL                   4                    [1, 1, 2, 2]
            CN                   4                    [1, 2, 2, 3]
            D?{                  5                    [1, 1, 1, 1, 4]
            D@s                  5                    [1, 1, 1, 2, 3]
            D@{                  5                    [1, 1, 2, 2, 4]
            DBg                  5                    [1, 1, 2, 2, 2]
            DBk                  5                    [1, 1, 2, 3, 3]
            DIk                  5                    [1, 2, 2, 2, 3]
            DK[                  5                    [1, 2, 2, 2, 3]
            D_K                  5                    [1, 1, 1, 1, 2]
            D`K                  5                    [1, 1, 2, 2, 2]
            E?Bw                 6                    [1, 1, 1, 1, 1, 5]
            E?Fg                 6                    [1, 1, 1, 1, 2, 4]
            E?N?                 6                    [1, 1, 1, 1, 2, 2]
            E?NG                 6                    [1, 1, 1, 1, 3, 3]
            E@FG                 6                    [1, 1, 1, 2, 2, 3]
            E@N?                 6                    [1, 1, 2, 2, 2, 2]
            E@Q?                 6                    [1, 1, 1, 1, 1, 1]
            E@QW                 6                    [1, 1, 1, 2, 2, 3]
            E@YO                 6                    [1, 1, 2, 2, 2, 2]
            E_?w                 6                    [1, 1, 1, 1, 1, 3]
            E_Cg                 6                    [1, 1, 1, 1, 2, 2]
            E_Cw                 6                    [1, 1, 1, 2, 2, 3]
            E_Ko                 6                    [1, 1, 2, 2, 2, 2]
            F??^?                7                    [1, 1, 1, 1, 1, 2, 3]
            F?LCG                7                    [1, 1, 1, 1, 2, 2, 2]
            FK??W                7                    [1, 1, 1, 1, 1, 1, 2]
            FK?GW                7                    [1, 1, 1, 1, 2, 2, 2]
            F_?@w                7                    [1, 1, 1, 1, 1, 1, 4]
            F_?Hg                7                    [1, 1, 1, 1, 1, 2, 3]
            F_?XO                7                    [1, 1, 1, 1, 2, 2, 2]
        """
        if graph_db is None:
            graph_db = GraphDatabase()
        if query_dict is not None:
            if query_dict['expression'][0] == 'degree_sequence':
                query_dict['expression'][3] = degseq_to_data(query_dict['expression'][3])
            elif query_dict['expression'][0] == 'induced_subgraphs':
                query_dict['expression'][3] = subgraphs_to_data(query_dict['expression'][3])
            SQLQuery.__init__(self, graph_db, query_dict)
        else:
            # construct a query from the given parameters
            SQLQuery.__init__(self, graph_db)

            # if display_cols is None:
            #    raise TypeError, 'Nonetype display_cols cannot retrieve data.'

            master_join = {}

            for key in kwds:
                # check validity
                if key not in valid_kwds:
                    raise KeyError('%s is not a valid key for this database.' % str(key))

                # designate a query_dict and reserve display_cols until end
                # (database.py currently concatenates them including repeats)
                qdict = {'display_cols': None}

                # set table name
                if key in graph_data:
                    qdict['table_name'] = 'graph_data'
                elif key in aut_grp:
                    qdict['table_name'] = 'aut_grp'
                elif key in degrees:
                    qdict['table_name'] = 'degrees'
                elif key in misc:
                    qdict['table_name'] = 'misc'
                elif key in spectrum:
                    qdict['table_name'] = 'spectrum'

                # set expression
                if not isinstance(kwds[key], list):
                    if key == 'induced_subgraphs':
                        s6 = (graph.Graph(kwds[key]).canonical_label()).graph6_string()
                        qdict['expression'] = [key, 'regexp', '.*%s.*' % s6]
                    else:
                        qdict['expression'] = [key, '=', kwds[key]]
                elif key == 'degree_sequence':
                    qdict['expression'] = [key, '=', degseq_to_data(kwds[key])]
                elif key != 'induced_subgraphs':
                    qdict['expression'] = [key] + kwds[key]

                # add key parameter to query
                join_dict = {qdict['table_name']: ('graph_id', 'graph_id')}
                if key == 'induced_subgraphs' and isinstance(kwds[key], list):
                    self.intersect(subgraphs_to_query(kwds[key], graph_db),
                                   'graph_data', join_dict, in_place=True)
                else:
                    self.intersect(SQLQuery(graph_db, qdict), 'graph_data',
                                   join_dict, in_place=True)

                # include search params (keys) in join clause
                # again, we exclude graph_data because it is the base table
                if qdict['table_name'] != 'graph_data':
                    master_join[qdict['table_name']] = ('graph_id', 'graph_id')

            # display columns from each table
            aut_grp_disp = ['aut_grp']
            degrees_disp = ['degrees']
            misc_disp = ['misc']
            spectrum_disp = ['spectrum']
            graph_data_disp = ['graph_data']

            disp_tables = [aut_grp_disp, degrees_disp, misc_disp, spectrum_disp]
            # graph_data intentionally left out because it is always called

            # organize display
            if display_cols is not None:
                for col in display_cols:
                    if col in graph_data:
                        graph_data_disp.append(col)
                    elif col in aut_grp:
                        aut_grp_disp.append(col)
                    elif col in degrees:
                        degrees_disp.append(col)
                    elif col in misc:
                        misc_disp.append(col)
                    elif col in spectrum:
                        spectrum_disp.append(col)

                # finish filling master join with display tables
                for tab in disp_tables:
                    if len(tab) > 1:
                        master_join[tab[0]] = ('graph_id', 'graph_id')

                # join clause for display tables
                join_str = 'FROM graph_data '
                for tab in master_join:
                    join_str += 'INNER JOIN %s ON graph_data.graph_id=%s.graph_id ' % (tab, tab)

                # construct sql syntax substring for display cols
                disp_list = ['SELECT graph_data.graph6, ']
                for col in graph_data_disp[1:]:
                    if col != 'graph6':
                        disp_list.append('graph_data.%s, ' % col)
                for col in aut_grp_disp[1:]:
                    disp_list.append('aut_grp.%s, ' % col)
                for col in degrees_disp[1:]:
                    disp_list.append('degrees.%s, ' % col)
                for col in misc_disp[1:]:
                    disp_list.append('misc.%s, ' % col)
                for col in spectrum_disp[1:]:
                    disp_list.append('spectrum.%s, ' % col)
                disp_list[-1] = disp_list[-1].rstrip(', ') + ' '
                disp_str = ''.join(disp_list)

                # substitute disp_str and join_str back into self's query string
                self.__query_string__ = re.sub('SELECT.*WHERE ',
                                               disp_str + join_str + 'WHERE ',
                                               self.__query_string__)
                self.__query_string__ += ' ORDER BY graph_data.graph6'

    def query_iterator(self):
        """
        Return an iterator over the results list of the :class:`~GraphQuery`.

        EXAMPLES::

            sage: Q = GraphQuery(display_cols=['graph6'], num_vertices=7, diameter=5)
            sage: for g in Q:
            ....:     print(g.graph6_string())
            F?`po
            F?gqg
            F@?]O
            F@OKg
            F@R@o
            FA_pW
            FEOhW
            FGC{o
            FIAHo
            sage: Q = GraphQuery(display_cols=['graph6'], num_vertices=7, diameter=5)
            sage: it = iter(Q)
            sage: while True:
            ....:     try: print(next(it).graph6_string())
            ....:     except StopIteration: break
            F?`po
            F?gqg
            F@?]O
            F@OKg
            F@R@o
            FA_pW
            FEOhW
            FGC{o
            FIAHo
        """
        return iter(self.get_graphs_list())

    __iter__ = query_iterator

    def show(self, max_field_size=20, with_picture=False):
        """
        Display the results of a query in table format.

        INPUT:

        - ``max_field_size`` -- integer (default: 20); width of fields in
          command prompt version

        - ``with_picture`` -- boolean (default: ``False``); whether or not to
          display results with a picture of the graph (available only in the
          notebook)

        EXAMPLES::

            sage: G = GraphDatabase()
            sage: Q = GraphQuery(G, display_cols=['graph6','num_vertices','aut_grp_size'], num_vertices=4, aut_grp_size=4)
            sage: Q.show()
            Graph6               Num Vertices         Aut Grp Size
            ------------------------------------------------------------
            C@                   4                    4
            C^                   4                    4

        ::

            sage: R = GraphQuery(G, display_cols=['graph6','num_vertices','degree_sequence'], num_vertices=4)
            sage: R.show()
            Graph6               Num Vertices         Degree Sequence
            ------------------------------------------------------------
            C?                   4                    [0, 0, 0, 0]
            C@                   4                    [0, 0, 1, 1]
            CB                   4                    [0, 1, 1, 2]
            CF                   4                    [1, 1, 1, 3]
            CJ                   4                    [0, 2, 2, 2]
            CK                   4                    [1, 1, 1, 1]
            CL                   4                    [1, 1, 2, 2]
            CN                   4                    [1, 2, 2, 3]
            C]                   4                    [2, 2, 2, 2]
            C^                   4                    [2, 2, 3, 3]
            C~                   4                    [3, 3, 3, 3]

        Show the pictures (in notebook mode only)::

            sage: S = GraphQuery(G, display_cols=['graph6','aut_grp_size'], num_vertices=4)
            sage: S.show(with_picture=True)
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot display plot on command line.

        Note that pictures can be turned off::

            sage: S.show(with_picture=False)
            Graph6               Aut Grp Size
            ----------------------------------------
            C?                   24
            C@                   4
            CB                   2
            CF                   6
            CJ                   6
            CK                   8
            CL                   2
            CN                   2
            C]                   8
            C^                   4
            C~                   24

        Show your own query (note that the output is not reformatted for
        generic queries)::

            sage: (GenericGraphQuery('select degree_sequence from degrees where max_degree=2 and min_degree >= 1', G)).show()
            degree_sequence
            --------------------
            211
            222
            2211
            2222
            21111
            22211
            22211
            22222
            221111
            221111
            222211
            222211
            222211
            222222
            222222
            2111111
            2221111
            2221111
            2221111
            2222211
            2222211
            2222211
            2222211
            2222222
            2222222
        """
        relabel = {}
        for col in valid_kwds:
            relabel[col] = ' '.join(word.capitalize() for word in col.split('_'))

        if re.search('SELECT .*degree_sequence.* FROM', self.__query_string__):
            format_cols = {'degree_sequence': data_to_degseq}
        else:
            format_cols = {}
        if with_picture:
            SQLQuery.show(self, max_field_size=max_field_size,
                          plot_cols={'graph6': graph6_to_plot},
                          format_cols=format_cols, id_col='graph6',
                          relabel_cols=relabel)
        else:
            SQLQuery.show(self, max_field_size=max_field_size,
                          format_cols=format_cols, relabel_cols=relabel,
                          id_col='graph6')

    def get_graphs_list(self):
        """
        Return a list of Sage Graph objects that satisfy the query.

        EXAMPLES::

            sage: Q = GraphQuery(display_cols=['graph6', 'num_vertices', 'degree_sequence'], num_edges=['<=', 5], min_degree=1)
            sage: L = Q.get_graphs_list()
            sage: L[0]
            Graph on 2 vertices
            sage: len(L)
            35
        """
        s = self.__query_string__
        re.sub('SELECT.*FROM ', 'SELECT graph6 FROM ', s)
        q = GenericGraphQuery(s, self.__database__, self.__param_tuple__)
        graph6_list = q.query_results()
        return [Graph(str(g[0])) for g in graph6_list]

    def number_of(self):
        """
        Return the number of graphs in the database that satisfy the query.

        EXAMPLES::

            sage: Q = GraphQuery(display_cols=['graph6', 'num_vertices', 'degree_sequence'] ,num_edges=['<=', 5], min_degree=1)
            sage: Q.number_of()
            35
        """
        # run graphs_list and return len
        s = self.__query_string__
        re.sub('SELECT.*FROM ', 'SELECT graph6 FROM ', s)
        q = GenericGraphQuery(s, self.__database__, self.__param_tuple__)
        return len(q.query_results())


class GraphDatabase(SQLDatabase):

    def __init__(self):
        """
        Graph Database.

        This class interfaces with the ``sqlite`` database ``graphs.db``. It is
        an immutable database that inherits from
        :class:`~sage.databases.sql_db.SQLDatabase` (see
        :mod:`sage.databases.sql_db`). The display functions and get_graphs_list
        create their own queries, but it is also possible to query the database
        by constructing either a :class:`~sage.databases.sql_db.SQLQuery`.

        The database contains all unlabeled graphs with 7 or fewer nodes. This
        class will also interface with the optional database package containing
        all unlabeled graphs with 8 or fewer nodes. The database consists of
        five tables. For a full table and column structure, call
        :func:`~graph_db_info`.

        The tables are associated by the unique primary key ``graph_id`` (int).

        To query this database, we create a :class:`~GraphQuery`. This can be
        done directly with the :meth:`~GraphDatabase.query` method or by
        initializing one of:

        - :class:`~GenericGraphQuery` -- allows direct entry of a query string
          and tuple of parameters. This is the route for more advanced users
          that are familiar with SQL

        - :class:`~GraphQuery` -- a wrapper of SQLQuery, a general
          database/query wrapper of SQLite for new users

        REFERENCES:

        - Data provided by Jason Grout (Brigham Young
          University). [Online] Available:
          http://artsci.drake.edu/grout/graphs/

        EXAMPLES::

            sage: G = GraphDatabase()
            sage: G.get_skeleton()
            {'aut_grp': {'aut_grp_size': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'edge_transitive': {'index': True,
               'primary_key': False,
               'sql': 'BOOLEAN',
               'unique': False},
              'graph_id': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'num_fixed_points': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'num_orbits': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'vertex_transitive': {'index': True,
               'primary_key': False,
               'sql': 'BOOLEAN',
               'unique': False}},
             'degrees': {'average_degree': {'index': True,
               'primary_key': False,
               'sql': 'REAL',
               'unique': False},
              'degree_sequence': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'degrees_sd': {'index': True,
               'primary_key': False,
               'sql': 'REAL',
               'unique': False},
              'graph_id': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'max_degree': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'min_degree': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'regular': {'index': True,
               'primary_key': False,
               'sql': 'BOOLEAN',
               'unique': False}},
             'graph_data': {'complement_graph6': {'index': True,
               'primary_key': False,
               'sql': 'TEXT',
               'unique': False},
              'eulerian': {'index': True,
               'primary_key': False,
               'sql': 'BOOLEAN',
               'unique': False},
              'graph6': {'index': True,
               'primary_key': False,
               'sql': 'TEXT',
               'unique': False},
              'graph_id': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': True},
              'lovasz_number': {'index': True,
               'primary_key': False,
               'sql': 'REAL',
               'unique': False},
              'num_cycles': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'num_edges': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'num_hamiltonian_cycles': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'num_vertices': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'perfect': {'index': True,
               'primary_key': False,
               'sql': 'BOOLEAN',
               'unique': False},
              'planar': {'index': True,
               'primary_key': False,
               'sql': 'BOOLEAN',
               'unique': False}},
             'misc': {'clique_number': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'diameter': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'edge_connectivity': {'index': True,
               'primary_key': False,
               'sql': 'BOOLEAN',
               'unique': False},
              'girth': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'graph_id': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'independence_number': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'induced_subgraphs': {'index': True,
               'primary_key': False,
               'sql': 'TEXT',
               'unique': False},
              'min_vertex_cover_size': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'num_components': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'num_cut_vertices': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'num_spanning_trees': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'radius': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'vertex_connectivity': {'index': True,
               'primary_key': False,
               'sql': 'BOOLEAN',
               'unique': False}},
             'spectrum': {'eigenvalues_sd': {'index': True,
               'primary_key': False,
               'sql': 'REAL',
               'unique': False},
              'energy': {'index': True,
               'primary_key': False,
               'sql': 'REAL',
               'unique': False},
              'graph_id': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'max_eigenvalue': {'index': True,
               'primary_key': False,
               'sql': 'REAL',
               'unique': False},
              'min_eigenvalue': {'index': True,
               'primary_key': False,
               'sql': 'REAL',
               'unique': False},
              'spectrum': {'index': False,
               'primary_key': False,
               'sql': 'TEXT',
               'unique': False}}}
        """
        SQLDatabase.__init__(self, dblocation)

    def _gen_interact_func(self, display, **kwds):
        """
        Return a function to interact with GraphQuery parameters and results.

        This is a helper method for the ``interactive_query`` method and should
        not be called directly.
        """
        function_name = '__temporary_interact_function'
        arg = ['%s=%s' % (word, kwds[word]) for word in kwds]
        boxes = ["%s=input_grid(1,2,['=',%s])" % (word, kwds[word]) for word in kwds]
        params = ['%s=%s[0]' % tuple(2 * [arg[i].split('=')[0]]) for i in range(len(arg))]

        s = 'def %s(%s):' % (function_name, ','.join(boxes))
        t = """
        print('<html><h2>Query Results:</h2></html>')
        GraphQuery(display_cols=%s,%s).show(with_picture=True)
        """ % tuple([display, ','.join(params)])
        s += '\t' + '\n\t'.join(t.split('\n')) + '\n'
        exec(s)
        return locals()[function_name]

    def query(self, query_dict=None, display_cols=None, **kwds):
        """
        Create a GraphQuery on this database.

        For full class details, type ``GraphQuery?``
        and press :kbd:`Shift` + :kbd:`Enter`.

        EXAMPLES::

            sage: D = GraphDatabase()
            sage: q = D.query(display_cols=['graph6', 'num_vertices', 'degree_sequence'], num_edges=['<=', 5])          # needs sage.symbolic
            sage: q.show()                                                              # needs sage.symbolic
            Graph6               Num Vertices         Degree Sequence
            ------------------------------------------------------------
            @                    1                    [0]
            A?                   2                    [0, 0]
            A_                   2                    [1, 1]
            B?                   3                    [0, 0, 0]
            BG                   3                    [0, 1, 1]
            BW                   3                    [1, 1, 2]
            Bw                   3                    [2, 2, 2]
            C?                   4                    [0, 0, 0, 0]
            C@                   4                    [0, 0, 1, 1]
            CB                   4                    [0, 1, 1, 2]
            CF                   4                    [1, 1, 1, 3]
            CJ                   4                    [0, 2, 2, 2]
            CK                   4                    [1, 1, 1, 1]
            CL                   4                    [1, 1, 2, 2]
            CN                   4                    [1, 2, 2, 3]
            C]                   4                    [2, 2, 2, 2]
            C^                   4                    [2, 2, 3, 3]
            D??                  5                    [0, 0, 0, 0, 0]
            D?C                  5                    [0, 0, 0, 1, 1]
            D?K                  5                    [0, 0, 1, 1, 2]
            D?[                  5                    [0, 1, 1, 1, 3]
            D?{                  5                    [1, 1, 1, 1, 4]
            D@K                  5                    [0, 0, 2, 2, 2]
            D@O                  5                    [0, 1, 1, 1, 1]
            D@S                  5                    [0, 1, 1, 2, 2]
            D@[                  5                    [0, 1, 2, 2, 3]
            D@s                  5                    [1, 1, 1, 2, 3]
            D@{                  5                    [1, 1, 2, 2, 4]
            DBW                  5                    [0, 2, 2, 2, 2]
            DB[                  5                    [0, 2, 2, 3, 3]
            DBg                  5                    [1, 1, 2, 2, 2]
            DBk                  5                    [1, 1, 2, 3, 3]
            DIk                  5                    [1, 2, 2, 2, 3]
            DK[                  5                    [1, 2, 2, 2, 3]
            DLo                  5                    [2, 2, 2, 2, 2]
            D_K                  5                    [1, 1, 1, 1, 2]
            D`K                  5                    [1, 1, 2, 2, 2]
            E???                 6                    [0, 0, 0, 0, 0, 0]
            E??G                 6                    [0, 0, 0, 0, 1, 1]
            E??W                 6                    [0, 0, 0, 1, 1, 2]
            E??w                 6                    [0, 0, 1, 1, 1, 3]
            E?@w                 6                    [0, 1, 1, 1, 1, 4]
            E?Bw                 6                    [1, 1, 1, 1, 1, 5]
            E?CW                 6                    [0, 0, 0, 2, 2, 2]
            E?C_                 6                    [0, 0, 1, 1, 1, 1]
            E?Cg                 6                    [0, 0, 1, 1, 2, 2]
            E?Cw                 6                    [0, 0, 1, 2, 2, 3]
            E?Dg                 6                    [0, 1, 1, 1, 2, 3]
            E?Dw                 6                    [0, 1, 1, 2, 2, 4]
            E?Fg                 6                    [1, 1, 1, 1, 2, 4]
            E?Ko                 6                    [0, 0, 2, 2, 2, 2]
            E?Kw                 6                    [0, 0, 2, 2, 3, 3]
            E?LO                 6                    [0, 1, 1, 2, 2, 2]
            E?LW                 6                    [0, 1, 1, 2, 3, 3]
            E?N?                 6                    [1, 1, 1, 1, 2, 2]
            E?NG                 6                    [1, 1, 1, 1, 3, 3]
            E@FG                 6                    [1, 1, 1, 2, 2, 3]
            E@HW                 6                    [0, 1, 2, 2, 2, 3]
            E@N?                 6                    [1, 1, 2, 2, 2, 2]
            E@Ow                 6                    [0, 1, 2, 2, 2, 3]
            E@Q?                 6                    [1, 1, 1, 1, 1, 1]
            E@QW                 6                    [1, 1, 1, 2, 2, 3]
            E@T_                 6                    [0, 2, 2, 2, 2, 2]
            E@YO                 6                    [1, 1, 2, 2, 2, 2]
            EG?W                 6                    [0, 1, 1, 1, 1, 2]
            EGCW                 6                    [0, 1, 1, 2, 2, 2]
            E_?w                 6                    [1, 1, 1, 1, 1, 3]
            E_Cg                 6                    [1, 1, 1, 1, 2, 2]
            E_Cw                 6                    [1, 1, 1, 2, 2, 3]
            E_Ko                 6                    [1, 1, 2, 2, 2, 2]
            F????                7                    [0, 0, 0, 0, 0, 0, 0]
            F???G                7                    [0, 0, 0, 0, 0, 1, 1]
            F???W                7                    [0, 0, 0, 0, 1, 1, 2]
            F???w                7                    [0, 0, 0, 1, 1, 1, 3]
            F??@w                7                    [0, 0, 1, 1, 1, 1, 4]
            F??Bw                7                    [0, 1, 1, 1, 1, 1, 5]
            F??GW                7                    [0, 0, 0, 0, 2, 2, 2]
            F??G_                7                    [0, 0, 0, 1, 1, 1, 1]
            F??Gg                7                    [0, 0, 0, 1, 1, 2, 2]
            F??Gw                7                    [0, 0, 0, 1, 2, 2, 3]
            F??Hg                7                    [0, 0, 1, 1, 1, 2, 3]
            F??Hw                7                    [0, 0, 1, 1, 2, 2, 4]
            F??Jg                7                    [0, 1, 1, 1, 1, 2, 4]
            F??Wo                7                    [0, 0, 0, 2, 2, 2, 2]
            F??Ww                7                    [0, 0, 0, 2, 2, 3, 3]
            F??XO                7                    [0, 0, 1, 1, 2, 2, 2]
            F??XW                7                    [0, 0, 1, 1, 2, 3, 3]
            F??Z?                7                    [0, 1, 1, 1, 1, 2, 2]
            F??ZG                7                    [0, 1, 1, 1, 1, 3, 3]
            F??^?                7                    [1, 1, 1, 1, 1, 2, 3]
            F?CJG                7                    [0, 1, 1, 1, 2, 2, 3]
            F?CPW                7                    [0, 0, 1, 2, 2, 2, 3]
            F?CZ?                7                    [0, 1, 1, 2, 2, 2, 2]
            F?C_w                7                    [0, 0, 1, 2, 2, 2, 3]
            F?Ca?                7                    [0, 1, 1, 1, 1, 1, 1]
            F?CaW                7                    [0, 1, 1, 1, 2, 2, 3]
            F?Ch_                7                    [0, 0, 2, 2, 2, 2, 2]
            F?CqO                7                    [0, 1, 1, 2, 2, 2, 2]
            F?LCG                7                    [1, 1, 1, 1, 2, 2, 2]
            F@??W                7                    [0, 0, 1, 1, 1, 1, 2]
            F@?GW                7                    [0, 0, 1, 1, 2, 2, 2]
            FG??w                7                    [0, 1, 1, 1, 1, 1, 3]
            FG?Gg                7                    [0, 1, 1, 1, 1, 2, 2]
            FG?Gw                7                    [0, 1, 1, 1, 2, 2, 3]
            FG?Wo                7                    [0, 1, 1, 2, 2, 2, 2]
            FK??W                7                    [1, 1, 1, 1, 1, 1, 2]
            FK?GW                7                    [1, 1, 1, 1, 2, 2, 2]
            F_?@w                7                    [1, 1, 1, 1, 1, 1, 4]
            F_?Hg                7                    [1, 1, 1, 1, 1, 2, 3]
            F_?XO                7                    [1, 1, 1, 1, 2, 2, 2]
        """
        return GraphQuery(self, query_dict, display_cols, **kwds)

    def interactive_query(self, display_cols, **kwds):
        """
        Generate an interact shell to query the database.

        .. WARNING::

            This is no longer implemented since the switch to Python3.

        This method generates an interact shell that allows the user
        to manipulate query parameters and see the updated results.

        .. TODO::

            This function could use improvement. Add full options of typical
            :class:`~GraphQuery` (i.e.: have it accept list input); and update
            options in interact to make it less annoying to put in operators.

        EXAMPLES::

            sage: D = GraphDatabase()
            sage: D.interactive_query(display_cols=['graph6', 'num_vertices', 'degree_sequence'], num_edges=5, max_degree=3)                                    # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: not available in Jupyter notebook
        """
        raise NotImplementedError('not available in Jupyter notebook')
        # print('<html><h1>Interactive Graph Query</h1></html>')
        # f = self._gen_interact_func(display=display_cols, **kwds)
        # interact(f)
