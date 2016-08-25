"""
Graph database

INFO:

This module implements classes (GraphDatabase, GraphQuery,
GenericGraphQuery) for interfacing with the sqlite database
graphs.db.

The GraphDatabase class interfaces with the sqlite database
graphs.db. It is an immutable database that inherits from
SQLDatabase (see sage.databases.database.py).

The database contains all unlabeled graphs with 7 or fewer nodes.
This class will also interface with the optional database package
containing all unlabeled graphs with 8 or fewer nodes. The
database(s) consists of five tables, and has the structure given by
the function graph_info. (For a full description including column
data types, create a GraphDatabase instance and call the method
get_skeleton).

AUTHORS:

- Emily A. Kirkman (2008-09-20): first version of interactive queries,
  cleaned up code and generalized many elements to
  sage.databases.database.py

- Emily A. Kirkman (2007-07-23): inherits GenericSQLDatabase, also
  added classes: GraphQuery and GenericGraphQuery

- Emily A. Kirkman (2007-05-11): initial sqlite version

- Emily A. Kirkman (2007-02-13): initial version (non-sqlite)

REFERENCES:

- Data provided by Jason Grout (Brigham Young University). [Online]
  Available: http://artsci.drake.edu/grout/graphs/
"""

################################################################################
#           Copyright (C) 2007 Emily A. Kirkman
#
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
################################################################################
from __future__ import print_function
from __future__ import absolute_import

from . import graph
import os
import re
from sage.rings.integer import Integer
from sqlite3 import dbapi2 as sqlite # if anyone would like to explain why dbapi2...
from sage.databases.sql_db import SQLDatabase, SQLQuery
from sage.env import SAGE_SHARE
from sage.graphs.graph import Graph
dblocation = os.path.join(SAGE_SHARE,'graphs','graphs.db')

def degseq_to_data(degree_sequence):
    """
    Takes a degree sequence list (of Integers) and converts to a sorted
    (max-min) integer data type, as used for faster access in the
    underlying database.

    EXAMPLE::

        sage: from sage.graphs.graph_database import degseq_to_data
        sage: degseq_to_data([2,2,3,1])
        3221
    """
    degree_sequence.sort()
    return sum(degree_sequence[i]*10**i for i in range(len(degree_sequence)))

def data_to_degseq(data, graph6=None):
    """
    Takes the database integer data type (one digit per vertex
    representing its degree, sorted high to low) and converts it to
    degree sequence list. The graph6 identifier is required for all
    graphs with no edges, so that the correct number of zeros will be
    returned.

    EXAMPLE::

        sage: from sage.graphs.graph_database import data_to_degseq
        sage: data_to_degseq(3221)
        [1, 2, 2, 3]
        sage: data_to_degseq(0,'D??')
        [0, 0, 0, 0, 0]
    """
    degseq = Integer(data).digits(10)
    if not degseq:
        # compute number of 0's in list from graph6 string
        from sage.graphs.generic_graph_pyx import length_and_string_from_graph6
        return length_and_string_from_graph6(str(graph6))[0]*[0]
    else:
        return degseq

def graph6_to_plot(graph6):
    """
    Constructs a graph from a graph6 string and returns a Graphics
    object with arguments preset for show function.

    EXAMPLE::

        sage: from sage.graphs.graph_database import graph6_to_plot
        sage: type(graph6_to_plot('D??'))
        <class 'sage.plot.graphics.Graphics'>
    """
    g = Graph(str(graph6))
    return g.plot(layout='circular',vertex_size=30,vertex_labels=False,graph_border=False)

def subgraphs_to_query(subgraphs, db):
    """
    Constructs and returns a GraphQuery object respecting the special
    input required for the induced_subgraphs parameter. This input can
    be an individual graph6 string (in which case it is evaluated
    without the use of this method) or a list of strings. In the latter
    case, the list should be of one of the following two formats: 1.
    ['one_of',String,...,String] Will search for graphs containing a
    subgraph isomorphic to any of the graph6 strings in the list. 2.
    ['all_of',String,...,String] Will search for graphs containing a
    subgraph isomorphic to each of the graph6 strings in the list.

    This is a helper method called by the GraphQuery constructor to
    handle this special format. This method should not be used on its
    own because it doesn't set any display columns in the query string,
    causing a failure to fetch the data when run.

    EXAMPLE::

        sage: from sage.graphs.graph_database import subgraphs_to_query
        sage: gd = GraphDatabase()
        sage: q = subgraphs_to_query(['all_of','A?','B?','C?'],gd)
        sage: q.get_query_string()
        'SELECT ,,,,,  FROM misc WHERE ( ( misc.induced_subgraphs regexp ? ) AND (
        misc.induced_subgraphs regexp ? ) ) AND ( misc.induced_subgraphs regexp ? )'
    """
    q = GraphQuery(graph_db=db,induced_subgraphs=subgraphs[1])
    if subgraphs[0] == 'all_of':
        for i in range(len(subgraphs))[2:]:
            q.intersect(GraphQuery(graph_db=db, induced_subgraphs=subgraphs[i]),in_place=True)
    elif subgraphs[0] == 'one_of':
        for i in range(len(subgraphs))[2:]:
            q.union(GraphQuery(graph_db=db, induced_subgraphs=subgraphs[i]),in_place=True)
    else:
        raise KeyError('Unable to initiate query:  Illegal input format for induced_subgraphs.')
    return q

# tables     columns                    input data type     sqlite data type
# -----------------------------------------------------------------------------
aut_grp =  ['aut_grp_size',             # Integer           INTEGER
            'num_orbits',               # Integer           INTEGER
            'num_fixed_points',         # Integer           INTEGER
            'vertex_transitive',        # bool              BOOLEAN
            'edge_transitive']          # bool              BOOLEAN
degrees =  ['degree_sequence',          # list              INTEGER (see degseq_to_data module function)
            'min_degree',               # Integer           INTEGER
            'max_degree',               # Integer           INTEGER
            'average_degree',           # Real              REAL
            'degrees_sd',               # Real              REAL
            'regular']                  # bool              BOOLEAN
misc =     ['vertex_connectivity',      # Integer           INTEGER
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
graph_data=['complement_graph6',        # String            STRING
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
    Returns a dictionary of allowed table and column names.

    INPUT:


    -  ``tablename`` - restricts the output to a single
       table


    EXAMPLE::

        sage: graph_db_info().keys()
        ['graph_data', 'degrees', 'spectrum', 'misc', 'aut_grp']

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
        A query for a GraphDatabase.

        INPUT:


        -  ``database`` - the GraphDatabase instance to query
           (if None then a new instance is created)

        -  ``query_string`` - a string representing the SQL
           query

        -  ``param_tuple`` - a tuple of strings - what to
           replace question marks in query_string with (optional, but a good
           idea)


        .. note::

           This query class is generally intended for developers and
           more advanced users. It allows you to execute any query,
           and so may be considered unsafe.

        See GraphDatabase class docstrings or enter::

            sage: G = GraphDatabase()
            sage: G.get_skeleton()
            {...

        to see the underlying structure of the database. Also see
        SQLQuery in sage.databases.database for more info and a
        tutorial.

        A piece of advice about '?' and param_tuple: It is generally
        considered safer to query with a '?' in place of each value
        parameter, and using a second argument (a tuple of strings) in a
        call to the sqlite database. Successful use of the param_tuple
        argument is exemplified::

            sage: G = GraphDatabase()
            sage: q = 'select graph_id,graph6,num_vertices,num_edges from graph_data where graph_id<=(?) and num_vertices=(?)'
            sage: param = (22,5)
            sage: Q = SQLQuery(G,q,param)
            sage: Q.show()
            graph_id             graph6               num_vertices         num_edges
            --------------------------------------------------------------------------------
            18                   D??                  5                    0
            19                   D?C                  5                    1
            20                   D?K                  5                    2
            21                   D@O                  5                    2
            22                   D?[                  5                    3
        """
        if database is None: database = GraphDatabase()
        if not isinstance(database, GraphDatabase):
            raise TypeError('%s is not a valid GraphDatabase'%database)
        SQLQuery.__init__(self,database,query_string,param_tuple)

class GraphQuery(GenericGraphQuery):

    def __init__(self, graph_db=None, query_dict=None, display_cols=None, **kwds):
        """
        A query for an instance of GraphDatabase. This class nicely wraps
        the SQLQuery class located in sage.databases.database.py to make
        the query constraints intuitive and with as many pre-definitions as
        possible. (i.e.: since it has to be a GraphDatabase, we already know
        the table structure and types; and since it is immutable, we can
        treat these as a guarantee).

        .. note::

           SQLQuery functions are available for GraphQuery. See
           sage.dataabases.database.py for more details.

        INPUT:


        -  ``graph_db`` - The GraphDatabase instance to apply
           the query to. (If None, then a new instance is created).

        -  ``query_dict`` - A dictionary specifying the query
           itself. Format is: 'table_name': 'tblname', 'display_cols':
           ['col1', 'col2'], 'expression':[col, operator, value] If not None,
           query_dict will take precedence over all other arguments.

        -  ``display_cols`` - A list of column names (strings)
           to display in the result when running or showing a query.

        -  ``kwds`` - The columns of the database are all
           keywords. For a database table/column structure dictionary, call
           graph_db_info. Keywords accept both single values and lists of
           length 2. The list allows the user to specify an expression other
           than equality. Valid expressions are strings, and for numeric
           values (i.e. Reals and Integers) are: '=','','','=','='. String
           values also accept 'regexp' as an expression argument. The only
           keyword exception to this format is induced_subgraphs, which
           accepts one of the following options: 1.
           ['one_of',String,...,String] Will search for graphs containing a
           subgraph isomorphic to any of the graph6 strings in the list. 2.
           ['all_of',String,...,String] Will search for graphs containing a
           subgraph isomorphic to each of the graph6 strings in the list.


        EXAMPLES::

            sage: Q = GraphQuery(display_cols=['graph6','num_vertices','degree_sequence'],num_edges=['<=',5],min_degree=1)
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
        if graph_db is None: graph_db = GraphDatabase()
        if query_dict is not None:
            if query_dict['expression'][0] == 'degree_sequence':
                query_dict['expression'][3] = degseq_to_data(query_dict['expression'][3])
            elif query_dict['expression'][0] == 'induced_subgraphs':
                query_dict['expression'][3] = subgraphs_to_data(query_dict['expression'][3])
            SQLQuery.__init__(self, graph_db, query_dict)
        else:
            # construct a query from the given parameters
            SQLQuery.__init__(self, graph_db)

            #if display_cols is None:
            #    raise TypeError, 'Nonetype display_cols cannot retrieve data.'

            master_join = {}

            for key in kwds:
                # check validity
                if not key in valid_kwds:
                    raise KeyError('%s is not a valid key for this database.'%str(key))

                # designate a query_dict
                qdict = {'display_cols': None}  # reserve display cols until end
                                                # (database.py currently concatenates
                                                # them including repeats)

                # set table name
                if key in graph_data: qdict['table_name'] = 'graph_data'
                elif key in aut_grp: qdict['table_name'] = 'aut_grp'
                elif key in degrees: qdict['table_name'] = 'degrees'
                elif key in misc: qdict['table_name'] = 'misc'
                elif key in spectrum: qdict['table_name'] = 'spectrum'

                # set expression
                if not isinstance(kwds[key],list):
                                        if key == 'induced_subgraphs':
                                                qdict['expression'] = [key, 'regexp', '.*%s.*'%(graph.Graph(kwds[key]).canonical_label()).graph6_string()]
                                        else:
                                                qdict['expression'] = [key, '=', kwds[key]]
                elif key == 'degree_sequence':
                    qdict['expression'] = [key, '=', degseq_to_data(kwds[key])]
                elif key != 'induced_subgraphs':
                    qdict['expression'] = [key] + kwds[key]

                # add key parameter to query
                join_dict = {qdict['table_name']: ('graph_id', 'graph_id')}
                if key == 'induced_subgraphs' and isinstance(kwds[key],list):
                    self.intersect(subgraphs_to_query(kwds[key], graph_db),
                                    'graph_data', join_dict, in_place=True)
                else:
                    self.intersect(SQLQuery(graph_db, qdict), 'graph_data', join_dict,in_place=True)

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
                                        if col in graph_data: graph_data_disp.append(col)
                                        elif col in aut_grp: aut_grp_disp.append(col)
                                        elif col in degrees: degrees_disp.append(col)
                                        elif col in misc: misc_disp.append(col)
                                        elif col in spectrum: spectrum_disp.append(col)

                                # finish filling master join with display tables
                                for tab in disp_tables:
                                        if len(tab) > 1:
                                                master_join[tab[0]] = ('graph_id', 'graph_id')

                                # join clause for display tables
                                join_str = 'FROM graph_data '
                                for tab in master_join:
                                        join_str += 'INNER JOIN %s ON graph_data.graph_id=%s.graph_id '%(tab, tab)

                                # construct sql syntax substring for display cols
                                disp_str = 'SELECT graph_data.graph6, '
                                for col in graph_data_disp[1:]:
                                        if col != 'graph6': disp_str += 'graph_data.%s, '%col
                                for col in aut_grp_disp[1:]: disp_str += 'aut_grp.%s, '%col
                                for col in degrees_disp[1:]: disp_str += 'degrees.%s, '%col
                                for col in misc_disp[1:]: disp_str += 'misc.%s, '%col
                                for col in spectrum_disp[1:]: disp_str += 'spectrum.%s, '%col
                                disp_str = disp_str.rstrip(', ') + ' '

                                # substitue disp_str and join_str back into self's query string
                                self.__query_string__ = re.sub('SELECT.*WHERE ', disp_str + join_str + \
                                                                                                'WHERE ', self.__query_string__)
                                self.__query_string__ += ' ORDER BY graph_data.graph6'

    def query_iterator(self):
        """
        Returns an iterator over the results list of the GraphQuery.

        EXAMPLE::

            sage: Q = GraphQuery(display_cols=['graph6'],num_vertices=7, diameter=5)
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
            sage: Q = GraphQuery(display_cols=['graph6'],num_vertices=7, diameter=5)
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
        Displays the results of a query in table format.

        INPUT:


        -  ``max_field_size`` - width of fields in command
           prompt version

        -  ``with_picture`` - whether or not to display
           results with a picture of the graph (available only in the
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

            sage: (GenericGraphQuery('select degree_sequence from degrees where max_degree=2 and min_degree >= 1',G)).show()
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
            relabel[col] = ' '.join([word.capitalize() for word in col.split('_')])

        if re.search('SELECT .*degree_sequence.* FROM',self.__query_string__):
            format_cols = {'degree_sequence': (lambda x,y: data_to_degseq(x,y))}
        else: format_cols = {}
        if with_picture:
            SQLQuery.show(self, max_field_size=max_field_size, \
                                plot_cols={'graph6': (lambda x: graph6_to_plot(x))}, \
                                format_cols=format_cols, id_col='graph6', \
                                relabel_cols=relabel)
        else:
            SQLQuery.show(self, max_field_size=max_field_size, \
                                format_cols=format_cols, relabel_cols=relabel, \
                                id_col='graph6')

    def get_graphs_list(self):
        """
        Returns a list of Sage Graph objects that satisfy the query.

        EXAMPLES::

            sage: Q = GraphQuery(display_cols=['graph6','num_vertices','degree_sequence'],num_edges=['<=',5],min_degree=1)
            sage: L = Q.get_graphs_list()
            sage: L[0]
            Graph on 2 vertices
            sage: len(L)
            35
        """
        from sage.graphs.graph_list import from_graph6

        s = self.__query_string__
        re.sub('SELECT.*FROM ', 'SELECT graph6 FROM ', s)
        q = GenericGraphQuery(s, self.__database__, self.__param_tuple__)
        graph6_list = q.query_results()
        return [Graph(str(g[0])) for g in graph6_list]

    def number_of(self):
        """
        Returns the number of graphs in the database that satisfy the
        query.

        EXAMPLES::

            sage: Q = GraphQuery(display_cols=['graph6','num_vertices','degree_sequence'],num_edges=['<=',5],min_degree=1)
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
        Graph Database

        INFO:

        This class interfaces with the sqlite database graphs.db. It is an
        immutable database that inherits from SQLDatabase (see
        sage.databases.database.py). The display functions and get_graphs_list
        create their own queries, but it is also possible to query the
        database by constructing either a SQLQuery.

        The database contains all unlabeled graphs with 7 or fewer nodes.
        This class will also interface with the optional database package
        containing all unlabeled graphs with 8 or fewer nodes. The database
        consists of five tables. For a full table and column structure,
        call graph_db_info.

        USE: The tables are associated by the unique primary key graph_id
        (int).

        To query this database, we create a GraphQuery. This can be done
        directly with the query method or by initializing one of 1.
        GenericGraphQuery - allows direct entry of a query string and tuple
        of parameters. This is the route for more advanced users that are
        familiar with SQL. 2. GraphQuery - is a wrapper of SQLQuery, a
        general database/query wrapper of SQLite for new users.

        REFERENCES:

        - Data provided by Jason Grout (Brigham Young
          University). [Online] Available:
          http://artsci.drake.edu/grout/graphs/

        EXAMPLE::

            sage: G = GraphDatabase()
            sage: G.get_skeleton()
            {u'aut_grp': {u'aut_grp_size': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'edge_transitive': {'index': True,
               'primary_key': False,
               'sql': u'BOOLEAN',
               'unique': False},
              u'graph_id': {'index': False,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'num_fixed_points': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'num_orbits': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'vertex_transitive': {'index': True,
               'primary_key': False,
               'sql': u'BOOLEAN',
               'unique': False}},
             u'degrees': {u'average_degree': {'index': True,
               'primary_key': False,
               'sql': u'REAL',
               'unique': False},
              u'degree_sequence': {'index': False,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'degrees_sd': {'index': True,
               'primary_key': False,
               'sql': u'REAL',
               'unique': False},
              u'graph_id': {'index': False,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'max_degree': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'min_degree': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'regular': {'index': True,
               'primary_key': False,
               'sql': u'BOOLEAN',
               'unique': False}},
             u'graph_data': {u'complement_graph6': {'index': True,
               'primary_key': False,
               'sql': u'TEXT',
               'unique': False},
              u'eulerian': {'index': True,
               'primary_key': False,
               'sql': u'BOOLEAN',
               'unique': False},
              u'graph6': {'index': True,
               'primary_key': False,
               'sql': u'TEXT',
               'unique': False},
              u'graph_id': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': True},
              u'lovasz_number': {'index': True,
               'primary_key': False,
               'sql': u'REAL',
               'unique': False},
              u'num_cycles': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'num_edges': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'num_hamiltonian_cycles': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'num_vertices': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'perfect': {'index': True,
               'primary_key': False,
               'sql': u'BOOLEAN',
               'unique': False},
              u'planar': {'index': True,
               'primary_key': False,
               'sql': u'BOOLEAN',
               'unique': False}},
             u'misc': {u'clique_number': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'diameter': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'edge_connectivity': {'index': True,
               'primary_key': False,
               'sql': u'BOOLEAN',
               'unique': False},
              u'girth': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'graph_id': {'index': False,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'independence_number': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'induced_subgraphs': {'index': True,
               'primary_key': False,
               'sql': u'TEXT',
               'unique': False},
              u'min_vertex_cover_size': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'num_components': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'num_cut_vertices': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'num_spanning_trees': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'radius': {'index': True,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'vertex_connectivity': {'index': True,
               'primary_key': False,
               'sql': u'BOOLEAN',
               'unique': False}},
             u'spectrum': {u'eigenvalues_sd': {'index': True,
               'primary_key': False,
               'sql': u'REAL',
               'unique': False},
              u'energy': {'index': True,
               'primary_key': False,
               'sql': u'REAL',
               'unique': False},
              u'graph_id': {'index': False,
               'primary_key': False,
               'sql': u'INTEGER',
               'unique': False},
              u'max_eigenvalue': {'index': True,
               'primary_key': False,
               'sql': u'REAL',
               'unique': False},
              u'min_eigenvalue': {'index': True,
               'primary_key': False,
               'sql': u'REAL',
               'unique': False},
              u'spectrum': {'index': False,
               'primary_key': False,
               'sql': u'TEXT',
               'unique': False}}}
        """
        SQLDatabase.__init__(self,dblocation)

    def _gen_interact_func(self, display, **kwds):
        """
        Generates and returns a function to interact with GraphQuery
        parameters and results. This is a helper method for the
        interactive_query method and should not be called directly.

        EXAMPLE::

            sage: D = GraphDatabase()
            sage: D.interactive_query(display_cols=['graph6','num_vertices','degree_sequence'],num_edges=['<=',5],max_degree=3)
            <html>...</html>

        ::

            sage: G = GraphDatabase()
            sage: f = G._gen_interact_func(display=['graph6'], num_vertices=3)
            sage: type(f)
            <type 'function'>
            sage: interact(f)
            <html>...
        """
        from sagenb.notebook.interact import input_grid
        arg=['%s=%s'%(word,kwds[word]) for word in kwds]
        boxes=["%s=input_grid(1,2,['=',%s])"%(word,kwds[word]) for word in kwds]
        params = ['%s=%s[0]'%tuple(2*[arg[i].split('=')[0]]) for i in range(len(arg))]

        s = 'def _(%s):'%','.join(boxes)
        t = """
        print('<html><h2>Query Results:</h2></html>')
        GraphQuery(display_cols=%s,%s).show(with_picture=True)
        """%tuple([display,','.join(params)])
        s += '\t'+'\n\t'.join(t.split('\n'))+'\n'
        exec(s)
        return _

    def query(self, query_dict=None, display_cols=None, **kwds):
        """
        Creates a GraphQuery on this database. For full class details, type
        GraphQuery? and press shift+enter.

        EXAMPLE::

            sage: D = GraphDatabase()
            sage: q = D.query(display_cols=['graph6','num_vertices','degree_sequence'],num_edges=['<=',5])
            sage: q.show()
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
        TODO: This function could use improvement. Add full options of
        typical GraphQuery (i.e.: have it accept list input); and update
        options in interact to make it less annoying to put in operators.

        Generates an interact shell (in the notebook only) that allows the
        user to manipulate query parameters and see the updated results.

        EXAMPLE::

            sage: D = GraphDatabase()
            sage: D.interactive_query(display_cols=['graph6','num_vertices','degree_sequence'],num_edges=5,max_degree=3)
            <html>...</html>
        """
        from sagenb.notebook.interact import interact
        print('<html><h1>Interactive Graph Query</h1></html>')
        f = self._gen_interact_func(display=display_cols,**kwds)
        interact(f)

