r"""
Relational (sqlite) Databases Module

INFO:

    This module implements classes (SQLDatabase and SQLQuery (pythonic
    implementation for the user with little or no knowledge of sqlite))
    that wrap the basic functionality of sqlite.

    Databases are constructed via a triple indexed dictionary called a
    skeleton. A skeleton should be constructed to fit the following format::

        | - skeleton -- a triple-indexed dictionary
        |     - outer key -- table name
        |         - inner key -- column name
        |             - inner inner key -- one of the following:
        |                 - ``primary_key`` -- boolean, whether column has been set as
                            primary key
        |                 - ``index`` -- boolean, whether column has been set as index
        |                 - ``unique`` -- boolean, whether column has been set as unique
        |                 - ``sql`` -- one of ``'TEXT'``, ``'BOOLEAN'``, ``'INTEGER'``,
                            ``'REAL'``, or other user defined type

    An example skeleton of a database with one table, that table with one
    column::

        {'table1':{'col1':{'primary_key':False, 'index':True, 'sql':'REAL'}}}

    SQLDatabases can also be constructed via the add, drop, and commit
    functions. The vacuum function is also useful for restoring hard disk space
    after a database has shrunk in size.

    A SQLQuery can be constructed by providing a query_dict, which is a
    dictionary with the following sample format::

        {'table_name': 'tblname', 'display_cols': ['col1', 'col2', 'col3'], 'expression':[col, operator, value]}

    Finally a SQLQuery also allows the user to directly input the query
    string for a database, and also supports the '?' syntax by allowing an
    argument for a tuple of parameters to query.

    For full details, please see the tutorial. sage.graphs.graph_database.py
    is an example of implementing a database class in Sage using this
    interface.

AUTHORS:

- R. Andrew Ohana (2011-07-16):  refactored and rewrote most of the code;
  merged the Generic classes into the non-Generic versions; changed the
  skeleton format to include a boolean indicating whether the column stores
  unique keys; changed the index names so as to avoid potential ambiguity

- Emily A. Kirkman (2008-09-20): added functionality to generate plots and
  reformat output in show

- Emily A. Kirkman and Robert L. Miller (2007-06-17): initial version
"""
# FUTURE TODOs (Ignore for now):
#    - order by clause in query strings
#    - delete from query containing joins
#    - add data by column
#    - wrap sqlalchemy
#    - create query interface (with interact)
#    - allow kwds arguments to SQLQuery (like GraphQuery)

# ****************************************************************************
#       Copyright (C) 2011 R. Andrew Ohana <andrew.ohana@gmail.com>
#       Copyright (C) 2007 Emily A. Kirkman
#                          Robert L. Miller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sqlite3 as sqlite
import os
import re
from sage.misc.temporary_file import tmp_filename
from sage.structure.sage_object import SageObject

sqlite_keywords = ['ABORT','ACTION','ADD','AFTER','ALL','ALTER','ANALYZE',
    'AND','AS','ASC','ATTACH','AUTOINCREMENT','BEFORE','BEGIN','BETWEEN','BY',
    'CASCADE','CASE','CAST','CHECK','COLLATE','COLUMN','COMMIT','CONFLICT',
    'CONSTRAINT','CREATE','CROSS','CURRENT_DATE','CURRENT_TIME',
    'CURRENT_TIMESTAMP','DATABASE','DEFAULT','DEFERRABLE','DEFERRED','DELETE',
    'DESC','DETACH','DISTINCT','DROP','EACH','ELSE','END','ESCAPE','EXCEPT',
    'EXCLUSIVE','EXISTS','EXPLAIN','FAIL','FOR','FOREIGN','FROM','FULL',
    'GLOB','GROUP','HAVING','IF','IGNORE','IMMEDIATE','IN','INDEX','INDEXED',
    'INITIALLY','INNER','INSERT','INSTEAD','INTERSECT','INTO','IS','ISNULL',
    'JOIN','KEY','LEFT','LIKE','LIMIT','MATCH','NATURAL','NO','NOT','NOTNULL',
    'NULL','OF','OFFSET','ON','OR','ORDER','OUTER','PLAN','PRAGMA','PRIMARY',
    'QUERY','RAISE','REFERENCES','REGEXP','REINDEX','RELEASE','RENAME',
    'REPLACE','RESTRICT','RIGHT','ROLLBACK','ROW','SAVEPOINT','SELECT','SET',
    'TABLE','TEMP','TEMPORARY','THEN','TO','TRANSACTION','TRIGGER','UNION',
    'UNIQUE','UPDATE','USING','VACUUM','VALUES','VIEW','VIRTUAL','WHEN',
    'WHERE']


def regexp(expr, item):
    """
    Function to define regular expressions in pysqlite.

    OUTPUT:

    - ``True`` if parameter ``item`` matches the regular expression
      parameter ``expr``
    - ``False`` otherwise (i.e.: no match)

    REFERENCES:

    - [Ha2005]_

    EXAMPLES::

        sage: from sage.databases.sql_db import regexp
        sage: regexp('.s.*','cs')
        True
        sage: regexp('.s.*','ccs')
        False
        sage: regexp('.s.*','cscccc')
        True
    """
    r = re.compile(expr)
    return r.match(item) is not None


def verify_type(type):
    """
    Verify that the specified ``type`` is one of the allowed strings.

    EXAMPLES::

        sage: from sage.databases.sql_db import verify_type
        sage: verify_type('INT')
        True
        sage: verify_type('int')
        True
        sage: verify_type('float')
        Traceback (most recent call last):
        ...
        TypeError: float is not a legal type.
    """
    types = ['INTEGER','INT','BOOLEAN','REAL','TEXT','BOOL','BLOB','NOTYPE']
    if type.upper() not in types:
        raise TypeError('%s is not a legal type.' % type)
    return True


def verify_column(col_dict):
    """
    Verify that ``col_dict`` is in proper format, and return a dict with
    default values filled in. Proper format::

        {'primary_key':False, 'index':False, 'unique': False, 'sql':'REAL'}

    EXAMPLES::

        sage: from sage.databases.sql_db import verify_column
        sage: col = {'sql':'BOOLEAN'}
        sage: verify_column(col)
        {'index': False, 'primary_key': False, 'sql': 'BOOLEAN', 'unique': False}
        sage: col = {'primary_key':True, 'sql':'INTEGER'}
        sage: verify_column(col)
        {'index': True, 'primary_key': True, 'sql': 'INTEGER', 'unique': True}
        sage: verify_column({})
        Traceback (most recent call last):
        ...
        ValueError: SQL type must be declared, e.g. {'sql':'REAL'}.
    """
    d = {}
    d['primary_key'] = col_dict.get('primary_key', False)
    d['index'] = col_dict.get('index', False) or d['primary_key']
    d['unique'] = col_dict.get('unique', False) or d['primary_key']
    if 'sql' not in col_dict:
        raise ValueError("SQL type must be declared, e.g. {'sql':'REAL'}.")
    if verify_type(col_dict['sql']):
        d['sql'] = col_dict['sql']
    return d


def verify_operator(operator):
    """
    Check that ``operator`` is one of the allowed strings.
    Legal operators include the following strings:

    - '='
    - '<='
    - '>='
    - '<'
    - '>'
    - '<>'
    - 'like'
    - 'regexp'
    - 'is null'
    - 'is not null'

    EXAMPLES::

        sage: from sage.databases.sql_db import verify_operator
        sage: verify_operator('<=')
        True
        sage: verify_operator('regexp')
        True
        sage: verify_operator('is null')
        True
        sage: verify_operator('not_an_operator')
        Traceback (most recent call last):
        ...
        TypeError: not_an_operator is not a legal operator.
    """
    binaries = ['=','<=','>=','like','<','>','<>','regexp']
    unaries = ['is null','is not null']
    if operator not in binaries and operator not in unaries:
        raise TypeError('%s is not a legal operator.' % operator)
    return True


def construct_skeleton(database):
    """
    Construct a database skeleton from the sql data.  The skeleton data
    structure is a triple indexed dictionary of the following format::

        | - skeleton -- a triple-indexed dictionary
        |   - outer key -- table name
        |     - inner key -- column name
        |       - inner inner key -- one of the following:
        |         - ``primary_key`` -- boolean, whether column has been set as
                    primary key
        |         - ``index`` -- boolean, whether column has been set as index
        |         - ``unique`` -- boolean, whether column has been set as unique
        |         - ``sql`` -- one of ``'TEXT'``, ``'BOOLEAN'``, ``'INTEGER'``,
                    ``'REAL'``, or other user defined type

    An example skeleton of a database with one table, that table with one
    column::

        {'table1':{'col1':{'primary_key':False, 'unique': True, 'index':True, 'sql':'REAL'}}}

    EXAMPLES::

        sage: G = SQLDatabase(GraphDatabase().__dblocation__, False)
        sage: from sage.databases.sql_db import construct_skeleton
        sage: sorted(construct_skeleton(G))
        ['aut_grp', 'degrees', 'graph_data', 'misc', 'spectrum']
    """
    skeleton = {}
    cur = database.__connection__.cursor()
    exe = cur.execute("SELECT name FROM sqlite_master WHERE TYPE='table'")
    for table in exe.fetchall():
        skeleton[table[0]] = {}
        exe1 = cur.execute("PRAGMA table_info(%s)" % table[0])
        for col in exe1.fetchall():
            if not col[2]:
                typ = 'NOTYPE'
            else:
                typ = col[2]
            skeleton[table[0]][col[1]] = {'sql':typ,
                'primary_key':(col[5] != 0), 'index':(col[5] != 0), 'unique': False}
        exe2 = cur.execute("PRAGMA index_list(%s)" % table[0])
        for col in exe2.fetchall():
            if col[1].find('sqlite') == -1:
                if os.path.basename(database.__dblocation__) == 'graphs.db':
                    name = col[1]
                else:
                    name = col[1][len(table[0])+3:]
                skeleton[table[0]][name]['index'] = True
                skeleton[table[0]][name]['unique'] = bool(col[2])
            else:
                name = cur.execute("PRAGMA index_info(%s)" % col[1])
                name = name.fetchone()[2]
                skeleton[table[0]][name]['unique'] = bool(col[2])
    return skeleton


p = 0


def _create_print_table(cur, col_titles, **kwds):
    r"""
    Create a nice printable table from the cursor given with the given
    column titles.

    KEYWORDS:

    - ``max_field_size`` -- how wide each field can be
    - ``format_cols`` -- dictionary that allows the user to specify the
      format of a column's output by supplying a function. The format of
      the dictionary is::

          {'column_name':(lambda x: format_function(x))}

    - ``plot_cols`` -- dictionary that allows the user to specify that a
      plot should be drawn by the object generated by a data slice. Note
      that plot kwds are permitted. The dictionary format is::

          {'column_name':((lambda x: plot_function(x)), **kwds)}

    - ``relabel_cols`` -- dictionary to specify a relabeling of column
      headers. The dictionary format is::

          {'table_name':{'old_col_name':'new_col_name'}}

    - ``id_col`` -- reference to a column that can be used as an object
      identifier for each row

    - ``html_table`` -- boolean that if ``True`` creates an html table instead of
      a print table. Always set to ``True`` in the notebook

    EXAMPLES::

        sage: DB = SQLDatabase()
        sage: DB.create_table('simon',{'a1':{'sql':'bool', 'primary_key':False}, 'b2':{'sql':'int'}})
        sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
        sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
        sage: from sage.databases.sql_db import _create_print_table
        sage: cur = r.__database__.__connection__.cursor()
        sage: exe = cur.execute(r.__query_string__, r.__param_tuple__)
        sage: _create_print_table(cur, [des[0] for des in cur.description])
        'a1                  \n--------------------\n0                  \n1                  \n1                  '
    """
    fcol_index = []
    pcol_index = []

    if 'format_cols' in kwds:
        fcol_map = []
        for col in kwds['format_cols']:
            fcol_map.append(kwds['format_cols'][col])
            fcol_index.append(col_titles.index(col))
    if 'plot_cols' in kwds:
        pcol_map = []
        for col in kwds['plot_cols']:
            pcol_map.append(kwds['plot_cols'][col])
            pcol_index.append(col_titles.index(col))

    max_field_size = kwds['max_field_size'] if 'max_field_size' in kwds \
                    else 20
    id_col_index = col_titles.index(kwds['id_col']) if 'id_col' in kwds \
                   else None

    if 'relabel_cols' in kwds:
        relabel_cols = kwds['relabel_cols']
        for i in range(len(col_titles)):
            try:
                col_titles[i] = relabel_cols[col_titles[i]]
            except KeyError:
                continue
    global p
    p = 0

    def row_str(row, html):
        f = 0
        global p
        cur_str = []
        for index in range(len(col_titles)):
            if index in pcol_index:
                if html:
                    plot = pcol_map[p % len(pcol_map)](row[index])
                    plot.save('%d.png' % p, figsize=[1,1])
                    field_val = '     <td bgcolor=white align=center> ' \
                        + '%s <br> <img src="cell://%d.png"> ' % (row[index],p) \
                        + '</td>\n'
                    p += 1
                else:
                    raise NotImplementedError('Cannot display plot on '
                                              'command line.')
            else:
                if index in fcol_index:
                    if id_col_index is None:
                        field_val = fcol_map[f](row[index])
                    else:
                        field_val = fcol_map[f](row[index], row[id_col_index])
                    f += 1
                else:
                    field_val = row[index]
                if html:
                    field_val = '     <td bgcolor=white align=center> ' \
                        + str(field_val) + ' </td>\n'
                else:
                    field_val = str(field_val).ljust(max_field_size)
            cur_str.append(field_val)
        return ' '.join(cur_str)

    if 'html_table' in kwds and kwds['html_table']:
        # Notebook Version
        ret = '<html><!--notruncate-->\n'
        ret += '  <table bgcolor=lightgrey cellpadding=0>\n'
        ret += '    <tr>\n      <td bgcolor=white align=center> '
        ret += (' </td>\n      <td bgcolor=white '
               + 'align=center> ').join(col_titles)
        ret += ' </td>\n    </tr>\n'
        ret += '\n'.join(['    <tr>\n ' + row_str(row, True) + '    </tr>'
               for row in cur])
        ret += '\n  </table>\n</html>'
    else:
        # Command Prompt Version
        ret = ' '.join([col.ljust(max_field_size) for col in col_titles])
        ret += '\n' + '-' * max_field_size * len(col_titles) + '\n'
        ret += '\n'.join([row_str(row, False) for row in cur])
    return ret


class SQLQuery(SageObject):
    def __init__(self, database, *args, **kwds):
        """
        A query for a SQLite database.

        INPUT:

        - ``database`` -- a SQLDatabase object
        - ``query_dict`` -- dictionary specifying the query itself. The
          format is::

              {'table_name':'tblname', 'display_cols':['col1', 'col2','col3'], 'expression': [col, operator, value]}

        NOTE:
            Every SQL type we are using is ultimately represented as a string,
            so if you wish to save actual strings to a database, you actually
            need to do something like: '"value"'.

        See the documentation of ``SQLDatabase`` for an introduction to using
        SQLite in Sage.

        EXAMPLES::

            sage: D = SQLDatabase()
            sage: D.create_table('simon',{'a1':{'sql':'bool', 'primary_key':False}, 'b2':{'sql':'int'}})
            sage: D.add_data('simon',[(0,0),(1,2),(2,4)])
            sage: r = SQLQuery(D, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 3]})
            sage: r.show()
            a1
            --------------------
            0
            1

        Test that :issue:`27562` is fixed::

            sage: D = SQLDatabase()
            sage: r = SQLQuery(D, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 3]})
            Traceback (most recent call last):
            ...
            ValueError: Database has no table simon
            sage: D.create_table('simon',{'a1':{'sql':'bool', 'primary_key':False}, 'b2':{'sql':'int'}})
            sage: D.create_table('simon',{'a1':{'sql':'bool', 'primary_key':False}, 'b2':{'sql':'int'}})
            Traceback (most recent call last):
            ...
            ValueError: Database already has a table named simon
            sage: SQLQuery(D, {'table_name':'simon', 'display_cols':['a1'], 'expression':['c1','>',2]})
            Traceback (most recent call last):
                ...
            ValueError: Table has no column c1
        """
        if not isinstance(database, SQLDatabase):
            raise TypeError('%s is not a valid SQLDatabase' % database)
        self.__database__ = database
        total_args = len(args) + len(kwds)
        if total_args == 0:
            self.__query_string__ = ''
            self.__param_tuple__ = tuple()
            self.__query_dict__ = {}
            return
        for x in args:
            if isinstance(x,dict):
                if 'query_dict' not in kwds:
                    kwds['query_dict'] = x
            elif isinstance(x, str):
                if 'query_string' not in kwds:
                    kwds['query_string'] = x
            elif isinstance(x, tuple):
                if 'param_tuple' not in kwds:
                    kwds['param_tuple'] = x
        if total_args > 2 or not ('query_dict' in kwds or
              'query_string' in kwds) or ('query_dict' in kwds and
              ('param_tuple' in kwds or 'query_string' in kwds)):
            raise ValueError('Query must be constructed with either a '
                + 'dictionary or a string and tuple')

        if 'query_dict' in kwds:
            query_dict = kwds['query_dict']
        else:
            self.__query_string__ = kwds['query_string']
            if 'param_tuple' in kwds:
                self.__param_tuple__ = tuple(str(x) for x in kwds['param_tuple'])
            else:
                self.__param_tuple__ = tuple()
            return
        if query_dict:
            skel = database.__skeleton__
            if query_dict['table_name'] not in skel:
                raise ValueError("Database has no table %s"
                    % query_dict['table_name'])
            table_name = query_dict['table_name']
            if query_dict['display_cols'] is not None:
                for column in query_dict['display_cols']:
                    if column not in skel[table_name]:
                        raise ValueError("Table has no column %s" % column)
            if query_dict['expression'][0] not in skel[table_name]:
                raise ValueError("Table has no column %s"
                    % query_dict['expression'][0])

            self.__query_dict__ = query_dict
            self.__param_tuple__ = (str(query_dict['expression'][2]),)
            verify_operator(query_dict['expression'][1])
            if query_dict['display_cols'] is None:
                self.__query_string__ = 'SELECT , FROM %s WHERE ' % table_name \
                    + '%s.%s ' % (table_name, query_dict['expression'][0]) \
                    + '%s ?' % query_dict['expression'][1]
            else:
                query_dict['display_cols'] = ['%s.%s' % (table_name, x)
                    for x in query_dict['display_cols']]
                self.__query_string__ = 'SELECT ' \
                    + ', '.join(query_dict['display_cols']) + ' FROM ' \
                    + '%s WHERE %s.' % (table_name, table_name) \
                    + '%s ' % query_dict['expression'][0] \
                    + '%s ?' % query_dict['expression'][1]
        else:
            self.__query_dict__ = {}
            self.__param_tuple__ = tuple()
            self.__query_string__ = ''

    def __repr__(self):
        """
        Override the print output to display useful info regarding the
        query.

        EXAMPLES::

            sage: G = GraphDatabase()
            sage: q = 'SELECT graph_id,graph6,num_vertices,num_edges FROM graph_data WHERE graph_id<=(?) AND num_vertices=(?)'
            sage: param = (22,5)
            sage: SQLQuery(G,q,param)
            Query for sql database: ...graphs.db
            Query string: SELECT graph_id,graph6,num_vertices,num_edges FROM
                graph_data WHERE graph_id<=(?) AND num_vertices=(?)
            Parameter tuple: ('22', '5')
            sage: r = 'SELECT graph6 FROM graph_data WHERE num_vertices<=3'
            sage: SQLQuery(G,r)
            Query for sql database: ...graphs.db
            Query string: SELECT graph6 FROM graph_data WHERE num_vertices<=3
        """
        if not self.__query_string__:
            return 'Empty query on %s.' % self.__database__.__dblocation__
        return "Query for sql database: %s" % self.__database__.__dblocation__ \
            + "\nQuery string: %s" % self.__query_string__ \
            + ("\nParameter tuple: %s" % str(self.__param_tuple__) if
            self.__param_tuple__ else "")

    def get_query_string(self):
        """
        Return a copy of the query string.

        EXAMPLES::

            sage: G = GraphDatabase()
            sage: q = 'SELECT graph_id,graph6,num_vertices,num_edges FROM graph_data WHERE graph_id<=(?) AND num_vertices=(?)'
            sage: param = (22,5)
            sage: SQLQuery(G,q,param).get_query_string()
            'SELECT graph_id,graph6,num_vertices,num_edges FROM graph_data
            WHERE graph_id<=(?) AND num_vertices=(?)'
            sage: r = 'SELECT graph6 FROM graph_data WHERE num_vertices<=3'
            sage: SQLQuery(G,r).get_query_string()
            'SELECT graph6 FROM graph_data WHERE num_vertices<=3'
        """
        from copy import copy
        return copy(self.__query_string__)

    def __iter__(self):
        """
        Return an iterator over the results of the query.

        EXAMPLES::

            sage: G = GraphDatabase()
            sage: q = 'SELECT graph_id,graph6 FROM graph_data WHERE num_vertices=(?)'
            sage: param = (5,)
            sage: Q = SQLQuery(G,q,param)
            sage: it = Q.__iter__()
            sage: next(it)
            (18, 'D??')
            sage: next(it)
            (19, 'D?C')
            sage: skip = [next(it) for _ in range(15)]
            sage: next(it)
            (35, 'DBk')
        """
        try:
            cur = self.__database__.__connection__.cursor()
            return cur.execute(self.__query_string__, self.__param_tuple__)
        except sqlite.OperationalError:
            raise RuntimeError('Failure to fetch query.')

    def query_results(self):
        """
        Run the query by executing the ``__query_string__``. Return the
        results of the query in a list.

        EXAMPLES::

            sage: G = GraphDatabase()
            sage: q = 'SELECT graph_id,graph6,num_vertices,num_edges FROM graph_data WHERE graph_id<=(?) AND num_vertices=(?)'
            sage: param = (22,5)
            sage: Q = SQLQuery(G,q,param)
            sage: Q.query_results()
            [(18, 'D??', 5, 0), (19, 'D?C', 5, 1), (20, 'D?K', 5, 2),
             (21, 'D@O', 5, 2), (22, 'D?[', 5, 3)]
            sage: R = SQLQuery(G,{'table_name':'graph_data', 'display_cols':['graph6'], 'expression':['num_vertices','=',4]})
            sage: R.query_results()
            [('C?',), ('C@',), ('CB',), ('CK',), ('CF',), ('CJ',),
             ('CL',), ('CN',), ('C]',), ('C^',), ('C~',)]
        """
        return list(self)

    def show(self, **kwds):
        """
        Display the result of the query in table format.

        KEYWORDS:

        - ``max_field_size`` -- how wide each field can be
        - ``format_cols`` -- dictionary that allows the user to specify the
          format of a column's output by supplying a function. The format of
          the dictionary is::

              {'column_name':(lambda x: format_function(x))}

        - ``plot_cols`` -- dictionary that allows the user to specify that a
          plot should be drawn by the object generated by a data slice. Note
          that plot kwds are permitted. The dictionary format is::

              {'column_name':((lambda x: plot_function(x)), **kwds)}

        - ``relabel_cols`` -- dictionary to specify a relabeling of column
          headers. The dictionary format is::

              {'table_name':{'old_col_name':'new_col_name'}}

        - ``id_col`` -- reference to a column that can be used as an object
          identifier for each row

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool', 'primary_key':False}, 'b2':{'sql':'int'}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: r.show()
            a1
            --------------------
            0
            1
            1
            sage: D = GraphDatabase()
            sage: from sage.graphs.graph_database import valid_kwds, data_to_degseq
            sage: relabel = {}
            sage: for col in valid_kwds:
            ....:     relabel[col] = ' '.join([word.capitalize() for word in col.split('_')])
            sage: q = GraphQuery(display_cols=['graph6','degree_sequence'], num_vertices=4)
            sage: SQLQuery.show(q, format_cols={'degree_sequence':(lambda x,y: data_to_degseq(x,y))}, relabel_cols=relabel, id_col='graph6')
            Graph6               Degree Sequence
            ----------------------------------------
            C?                   [0, 0, 0, 0]
            C@                   [0, 0, 1, 1]
            CB                   [0, 1, 1, 2]
            CF                   [1, 1, 1, 3]
            CJ                   [0, 2, 2, 2]
            CK                   [1, 1, 1, 1]
            CL                   [1, 1, 2, 2]
            CN                   [1, 2, 2, 3]
            C]                   [2, 2, 2, 2]
            C^                   [2, 2, 3, 3]
            C~                   [3, 3, 3, 3]
        """
        if not self.__query_string__:
            return self.__database__.show()

        try:
            cur = self.__database__.__connection__.cursor()
            cur.execute(self.__query_string__, self.__param_tuple__)
        except Exception:
            raise RuntimeError('Failure to fetch query.')

        print(_create_print_table(cur, [des[0] for des in cur.description],
                                  **kwds))

    def __copy__(self):
        """
        Return a copy of itself.

        EXAMPLES::

            sage: G = GraphDatabase()
            sage: Q = SQLQuery(G, {'table_name':'graph_data', 'display_cols':['graph_id','graph6','num_vertices'], 'expression':['num_edges','<',3]})
            sage: R = copy(Q)
            sage: R.__query_string__ = ''
            sage: Q.__query_string__ == ''
            False
        """
        d = SQLQuery(self.__database__)
        d.__query_dict__ = self.__query_dict__
        d.__query_string__ = self.__query_string__
        d.__param_tuple__ = self.__param_tuple__
        return d

    def intersect(self, other, join_table=None, join_dict=None,
                  in_place=False):
        """
        Return a new ``SQLQuery`` that is the intersection of ``self`` and
        ``other``. ``join_table`` and ``join_dict`` can be ``None`` iff the
        two queries only search one table in the database. All display columns
        will be concatenated in order: ``self`` display cols + other display cols.

        INPUT:

        - ``other`` -- the ``SQLQuery`` to intersect with
        - ``join_table`` -- base table to join on (This table should have at
          least one column in each table to join on).
        - ``join_dict`` -- dictionary that represents the join structure for
          the new query. (Must include a mapping for all tables, including
          those previously joined in either query). Structure is given by::

              {'join_table1':('corr_base_col1', 'col1'), 'join_table2':('corr_base_col2', 'col2')}

          where ``join_table1`` is to be joined with ``join_table`` on
          ``join_table.corr_base_col1 = join_table1.col1``

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.create_table('lucy',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_data('simon', [(0,5),(1,4)])
            sage: DB.add_data('lucy', [(1,1),(1,4)])
            sage: q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['b2'], 'expression':['a1','=',1]})
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: s = q.intersect(r, 'simon', {'lucy':('a1','a1')})
            sage: s.get_query_string()
            'SELECT lucy.b2,simon.a1 FROM simon INNER JOIN lucy ON
            simon.a1=lucy.a1 WHERE ( lucy.a1 = ? ) AND ( simon.b2 <= ? )'
            sage: s.query_results()
            [(1, 1), (4, 1)]
            sage: s = q.intersect(r)
            Traceback (most recent call last):
            ...
            ValueError: Input queries query different tables but join
                parameters are NoneType
            sage: s.__query_string__ == q.__query_string__
            False
            sage: q.intersect(r, 'simon', {'lucy':('a1','a1')}, True)
            sage: q.__query_string__ == s.__query_string__
            True
        """
        if self.__query_dict__ is None or other.__query_dict__ is None:
            raise RuntimeError('Queries must be constructed using a '
                + 'dictionary in order to be intersected.')
        if self.__database__ != other.__database__:
            raise TypeError('Queries %s and %s must be ' % (self, other)
                + 'attached to the same database.')

        if in_place:
            if not self.__query_string__:
                self.__query_string__ = other.__query_string__
                self.__param_tuple__ = other.__param_tuple__
            elif not other.__query_string__:
                return
            else:
                self._merge_queries(other, self, join_table, join_dict, 'AND')
        else:
            from copy import copy
            if not self.__query_string__:
                return copy(other)
            if not other.__query_string__:
                return copy(self)
            return self._merge_queries(other, copy(self), join_table,
                join_dict, 'AND')

    def _merge_queries(self, other, ret, join_table, join_dict, operator):
        """
        The query ``ret`` is set to a new ``SQLQuery`` that is combines self
        and other through ``operator``. The other arguments are the same as
        ``intersect`` and ``union``.

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.create_table('lucy',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_data('simon', [(0,5),(1,4)])
            sage: DB.add_data('lucy', [(1,1),(1,4)])
            sage: q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['b2'], 'expression':['a1','=',1]})
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: s = q._merge_queries(r, copy(q), 'simon', {'lucy':('a1','a1')}, 'OR')
            sage: s.get_query_string()
            'SELECT lucy.b2,simon.a1 FROM simon INNER JOIN lucy ON
            simon.a1=lucy.a1 WHERE ( lucy.a1 = ? ) OR ( simon.b2 <= ? )'
            sage: s.query_results()
            [(1, 1), (4, 1)]
        """
        if join_table is None or join_dict is None:
            pattern = ' JOIN '
            if re.search(pattern, self.__query_string__) \
              or re.search(pattern, other.__query_string__):
                raise TypeError('Input queries have joins but join '
                    + 'parameters are NoneType')
            s = ((self.__query_string__).upper()).split('FROM ')
            o = ((other.__query_string__).upper()).split('FROM ')
            s = s[1].split(' WHERE ')
            o = o[1].split(' WHERE ')
            if s[0] != o[0]:
                raise ValueError('Input queries query different tables but '
                    + 'join parameters are NoneType')

        # inner join clause
        if join_dict is not None:
            joins = join_table
            for table in join_dict:
                joins += ' INNER JOIN %s ON %s.' % (table, join_table) \
                    + '%s=%s.' % (join_dict[table][0], table) \
                    + '%s ' % join_dict[table][1]
            ret.__query_string__ = re.sub(' FROM .* WHERE ', ' FROM ' + joins
                + 'WHERE ', self.__query_string__)

        # concatenate display cols
        disp1 = ret.__query_string__.split(' FROM')
        disp2 = other.__query_string__.split(' FROM')
        disp1.insert(1, ',%s FROM' % disp2[0].split('SELECT ')[1])
        new_query = ''.join(disp1)

        # concatenate where clause
        new_query = re.sub(' WHERE ', ' WHERE ( ',
                           new_query)
        new_query += re.sub('^.* WHERE ', f' ) {operator} ( ',
                            other.__query_string__)
        ret.__query_string__ = new_query + ' )'

        ret.__param_tuple__ = self.__param_tuple__ + other.__param_tuple__

        return ret

    def union(self, other, join_table=None, join_dict=None, in_place=False):
        """
        Return a new ``SQLQuery`` that is the union of ``self`` and ``other``.
        ``join_table`` and ``join_dict`` can be ``None`` iff the two queries
        only search one table in the database. All display columns will be
        concatenated in order: ``self`` display cols + other display cols.

        INPUT:

        - ``other`` -- the ``SQLQuery`` to union with
        - ``join_table`` -- base table to join on (This table should have at
          least one column in each table to join on).
        - ``join_dict`` -- dictionary that represents the join structure for
          the new query. (Must include a mapping for all tables, including
          those previously joined in either query). Structure is given by::

              {'join_table1':('corr_base_col1', 'col1'), 'join_table2':('corr_base_col2', 'col2')}

          where ``join_table1`` is to be joined with ``join_table`` on
          ``join_table.corr_base_col1=join_table1.col1``

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.create_table('lucy',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_data('simon', [(0,5),(1,4)])
            sage: DB.add_data('lucy', [(1,1),(1,4)])
            sage: q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['b2'], 'expression':['a1','=',1]})
            sage: r = SQLQuery(DB, {'table_name':'simon', 'display_cols':['a1'], 'expression':['b2','<=', 6]})
            sage: s = q.union(r, 'simon', {'lucy':('a1','a1')})
            sage: s.get_query_string()
            'SELECT lucy.b2,simon.a1 FROM simon INNER JOIN lucy ON
            simon.a1=lucy.a1 WHERE ( lucy.a1 = ? ) OR ( simon.b2 <= ? )'
            sage: s.query_results()
            [(1, 1), (4, 1)]
        """
        if self.__query_dict__ is None or other.__query_dict__ is None:
            raise RuntimeError('Queries must be constructed using a '
                + 'dictionary in order to be unioned.')
        if self.__database__ != other.__database__:
            raise TypeError('Queries %s and %s must be ' % (self, other)
                + 'attached to the same database.')

        if in_place:
            if self.__query_string__ and other.__query_string__:
                self._merge_queries(other, self, join_table, join_dict, 'OR')
        else:
            from copy import copy
            if not self.__query_string__:
                return copy(self)
            if not other.__query_string__:
                return copy(other)
            return self._merge_queries(other, copy(self), join_table,
                join_dict, 'OR')


class SQLDatabase(SageObject):
    def __init__(self, filename=None, read_only=None, skeleton=None):
        r"""
        A SQL Database object corresponding to a database file.

        INPUT:

        - ``filename`` -- string
        - ``skeleton`` -- a triple-indexed dictionary::

            | - outer key -- table name
            |   - inner key -- column name
            |     - inner inner key -- one of the following:
            |       - ``primary_key`` -- boolean, whether column has been set
                      as primary key
            |       - ``index`` -- boolean, whether column has been set as
                      index
            |       - ``unique`` -- boolean, whether column has been set as
                      unique
            |       - ``sql`` -- one of ``'TEXT'``, ``'BOOLEAN'``,
                      ``'INTEGER'``, ``'REAL'``, or other user defined type

        TUTORIAL:

        The ``SQLDatabase`` class is for interactively building databases
        intended for queries. This may sound redundant, but it is important. If
        you want a database intended for quick lookup of entries in very large
        tables, much like a hash table (such as a Python dictionary), a
        ``SQLDatabase`` may not be what you are looking for. The strength of
        ``SQLDatabases`` is in queries, searches through the database with
        complicated criteria.

        For example, we create a new database for storing isomorphism classes
        of simple graphs::

            sage: D = SQLDatabase()

        In order to generate representatives for the classes, we will import a
        function which generates all labeled graphs (noting that this is not
        the optimal way)::

            sage: from sage.groups.perm_gps.partn_ref.refinement_graphs import all_labeled_graphs

        We will need a table in the database in which to store the graphs, and
        we specify its structure with a Python dictionary, each of whose keys
        is the name of a column::

            sage: from collections import OrderedDict
            sage: table_skeleton = OrderedDict([
            ....: ('graph6',{'sql':'TEXT', 'index':True, 'primary_key':True}),
            ....: ('vertices', {'sql':'INTEGER'}),
            ....: ('edges', {'sql':'INTEGER'})
            ....: ])

        Then we create the table::

            sage: D.create_table('simon', table_skeleton)
            sage: D.show('simon')
            graph6               vertices             edges
            ------------------------------------------------------------

        Now that we have the table, we will begin to populate the table with
        rows. First, add the graph on zero vertices.::

            sage: G = Graph()
            sage: D.add_row('simon',(G.graph6_string(), 0, 0))
            sage: D.show('simon')
            graph6               vertices             edges
            ------------------------------------------------------------
            ?                    0                    0

        Next, add the graph on one vertex.::

            sage: G.add_vertex()
            0
            sage: D.add_row('simon',(G.graph6_string(), 1, 0))
            sage: D.show('simon')
            graph6               vertices             edges
            ------------------------------------------------------------
            ?                    0                    0
            @                    1                    0

        Say we want a database of graphs on four or less vertices::

            sage: labels = {}
            sage: for i in range(2, 5):
            ....:     labels[i] = []
            ....:     for g in all_labeled_graphs(i):
            ....:         g = g.canonical_label(algorithm='sage')
            ....:         if g not in labels[i]:
            ....:             labels[i].append(g)
            ....:             D.add_row('simon', (g.graph6_string(), g.order(), g.size()))
            sage: D.show('simon')
            graph6               vertices             edges
            ------------------------------------------------------------
            ?                    0                    0
            @                    1                    0
            A?                   2                    0
            A_                   2                    1
            B?                   3                    0
            BG                   3                    1
            BW                   3                    2
            Bw                   3                    3
            C?                   4                    0
            C@                   4                    1
            CB                   4                    2
            CF                   4                    3
            CJ                   4                    3
            CK                   4                    2
            CL                   4                    3
            CN                   4                    4
            C]                   4                    4
            C^                   4                    5
            C~                   4                    6

        We can then query the database -- let's ask for all the graphs on four
        vertices with three edges. We do so by creating two queries and asking
        for rows that satisfy them both::

            sage: Q = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['vertices','=',4]})
            sage: Q2 = SQLQuery(D, {'table_name':'simon', 'display_cols':['graph6'], 'expression':['edges','=',3]})
            sage: Q = Q.intersect(Q2)
            sage: len(Q.query_results())
            3
            sage: Q.query_results() # random
            [('CF', 'CF'), ('CJ', 'CJ'), ('CL', 'CL')]

        NOTE: The values of ``display_cols`` are always concatenated in
        intersections and unions.

        Of course, we can save the database to file. Here we use a
        temporary directory that we clean up afterwards::

            sage: import tempfile
            sage: d = tempfile.TemporaryDirectory()
            sage: dbpath = os.path.join(d.name, 'simon.db')
            sage: D.save(dbpath)

        Now the database's hard link is to this file, and not the temporary db
        file. For example, let's say we open the same file with another class
        instance. We can load the file as an immutable database::

            sage: E = SQLDatabase(dbpath)
            sage: E.show('simon')
            graph6               vertices             edges
            ------------------------------------------------------------
            ?                    0                    0
            @                    1                    0
            A?                   2                    0
            A_                   2                    1
            B?                   3                    0
            BG                   3                    1
            BW                   3                    2
            Bw                   3                    3
            C?                   4                    0
            C@                   4                    1
            CB                   4                    2
            CF                   4                    3
            CJ                   4                    3
            CK                   4                    2
            CL                   4                    3
            CN                   4                    4
            C]                   4                    4
            C^                   4                    5
            C~                   4                    6
            sage: E.drop_table('simon')
            Traceback (most recent call last):
            ...
            RuntimeError: Cannot drop tables from a read only database.

        Call ``cleanup()`` on the temporary directory to, well, clean it up::

            sage: d.cleanup()
        """
        if filename is None:
            if read_only is None:
                read_only = False
            filename = tmp_filename() + '.db'
        elif (filename[-3:] != '.db'):
            raise ValueError('Please enter a valid database path (file name '
                + '%s does not end in .db).' % filename)
        if read_only is None:
            read_only = True

        self.__read_only__ = read_only
        self.ignore_warnings = False
        self.__dblocation__ = filename
        self.__connection__ = sqlite.connect(self.__dblocation__,
            check_same_thread=False)
        # this is to avoid the multiple thread problem with dsage:
        # pysqlite does not trust multiple threads for the same connection
        self.__connection__.create_function("regexp", 2, regexp)

        # construct skeleton (from provided database)
        self.__skeleton__ = construct_skeleton(self)

        # add bones from new skeleton to database,
        # without changing existing structure
        if skeleton is not None and not read_only:
            for table in skeleton:
                if table not in self.__skeleton__:
                    self.create_table(table, skeleton[table])
                else:
                    for column in skeleton[table]:
                        if column not in self.__skeleton__[table]:
                            self.add_column(table, column,
                                            skeleton[table][column])
                        else:
                            print('Column attributes were ignored for '
                                'table {}, column {} -- column is '
                                'already in table.'.format(table, column))
        elif skeleton is not None:
            raise RuntimeError('Cannot update skeleton of a read only '
                + 'database.')

    def __repr__(self):
        """
        Override the print output to display useful info regarding the
        database.

        EXAMPLES::

            sage: import tempfile
            sage: with tempfile.TemporaryDirectory() as d:
            ....:     dbpath = os.path.join(d, "test.db")
            ....:     SD = SQLDatabase(dbpath, False)
            ....:     SD.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}})
            ....:     print(SD)
            table simon:
                column n: index: True; primary_key: False; sql: INTEGER;
                    unique: False;
        """
        s = ''
        for table in self.__skeleton__:
            s += 'table ' + table + ':\n'
            for column in self.__skeleton__[table]:
                s += '    column ' + column + ': '
                for data in sorted(self.__skeleton__[table][column]):
                    s += data + ': ' \
                        + str(self.__skeleton__[table][column][data]) + '; '
                s += '\n'
        return s

    def __copy__(self):
        """
        Return an instance of ``SQLDatabase`` that points to a copy database,
        and allows modification.

        EXAMPLES::

            sage: from collections import OrderedDict
            sage: DB = SQLDatabase()
            sage: DB.create_table('lucy', OrderedDict([
            ....: ('id', {'sql':'INTEGER', 'primary_key':True, 'index':True}),
            ....: ('a1', {'sql':'bool'}),
            ....: ('b2', {'sql':'int', 'primary_key':False})
            ....: ]))
            sage: DB.add_rows('lucy', [(0,1,1),(1,1,4),(2,0,7),(3,1,384), (4,1,978932)],['id','a1','b2'])
            sage: d = copy(DB)
            sage: d == DB
            False
            sage: d.show('lucy')
            id                   a1                   b2
            ------------------------------------------------------------
            0                    1                    1
            1                    1                    4
            2                    0                    7
            3                    1                    384
            4                    1                    978932
            sage: DB.show('lucy')
            id                   a1                   b2
            ------------------------------------------------------------
            0                    1                    1
            1                    1                    4
            2                    0                    7
            3                    1                    384
            4                    1                    978932
        """
        # copy .db file
        new_loc = tmp_filename() + '.db'
        if not self.__read_only__:
            self.commit()
        os.system('cp ' + self.__dblocation__ + ' ' + new_loc)
        D = SQLDatabase(filename=new_loc, read_only=False)
        return D

    def save(self, filename):
        """
        Save the database to the specified location.

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}})
            sage: for n in range(20): MonicPolys.add_row('simon', (n,))
            sage: import tempfile
            sage: with tempfile.TemporaryDirectory() as d:
            ....:     dbpath = os.path.join(d, "sage.db")
            ....:     MonicPolys.save(dbpath)
            ....:     N = SQLDatabase(dbpath)
            ....:     N.show('simon')
            n
            --------------------
            0
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
            13
            14
            15
            16
            17
            18
            19
        """
        if not self.__read_only__:
            self.commit()
        os.system('cp ' + self.__dblocation__ + ' ' + filename)

    def get_skeleton(self, check=False):
        """
        Return a dictionary representing the hierarchical structure of the
        database, in the following format::

            | - skeleton -- a triple-indexed dictionary
            |   - outer key -- table name
            |     - inner key -- column name
            |       - inner inner key -- one of the following:
            |         - ``primary_key`` -- boolean, whether column has been set as
                        primary key
            |         - ``index`` -- boolean, whether column has been set as index
            |         - ``unique`` -- boolean, whether column has been set as unique
            |         - ``sql`` -- one of ``'TEXT'``, ``'BOOLEAN'``, ``'INTEGER'``,
                        ``'REAL'``, or other user defined type

        For example::

            {'table1':{'col1':{'primary_key':False, 'index':True, 'unique': False,'sql':'REAL'}}}

        INPUT:

        - ``check`` -- if ``True``, checks to make sure the database's actual
          structure matches the skeleton on record

        EXAMPLES::

            sage: GDB = GraphDatabase()
            sage: GDB.get_skeleton()             # slightly random output
            {'aut_grp': {'aut_grp_size': {'index': True,
                                           'unique': False,
                                           'primary_key': False,
                                           'sql': 'INTEGER'},
                         ...
                         'num_vertices': {'index': True,
                                          'unique': False,
                                          'primary_key': False,
                                          'sql': 'INTEGER'}}}
        """
        if check:
            d = construct_skeleton(self)
            if d == self.__skeleton__:
                return d
            raise RuntimeError("skeleton structure is out of whack")
        return self.__skeleton__

    def query(self, *args, **kwds):
        """
        Create a ``SQLQuery`` on this database.

        For full class details,
        type ``SQLQuery?`` and press :kbd:`Shift` + :kbd:`Enter`.

        EXAMPLES::

            sage: D = SQLDatabase()
            sage: D.create_table('simon', {'wolf':{'sql':'BOOLEAN'}, 'tag':{'sql':'INTEGER'}})
            sage: q = D.query({'table_name':'simon', 'display_cols':['tag'], 'expression':['wolf','=',1]})
            sage: q.get_query_string()
            'SELECT simon.tag FROM simon WHERE simon.wolf = ?'
            sage: q.__param_tuple__
            ('1',)
            sage: q = D.query(query_string='SELECT tag FROM simon WHERE wolf=?',param_tuple=(1,))
            sage: q.get_query_string()
            'SELECT tag FROM simon WHERE wolf=?'
            sage: q.__param_tuple__
            ('1',)
        """
        return SQLQuery(self, *args, **kwds)

    __call__ = query

    def show(self, table_name, **kwds):
        """
        Show an entire table from the database.

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: DB.show('simon')
            a1                   b2
            ----------------------------------------
            0                    0
            1                    1
            1                    2
        """
        try:
            cur = self.__connection__.cursor()
            cur.execute('SELECT * FROM ' + table_name)
        except Exception:
            raise RuntimeError('Failure to fetch data.')
        print(_create_print_table(cur, [des[0] for des in cur.description],
                                  **kwds))

    def get_cursor(self, ignore_warning=None):
        """
        Return a pysqlite cursor for the database connection.

        A cursor is an input from which you can execute sqlite commands on the
        database.

        Recommended for more advanced users only.

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_row('simon',(0,1))
            sage: DB.add_rows('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_rows('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(0, 1), (0, 0), (1, 1), (1, 2), (0, 0), (0, 4), (1, 5)]
        """
        if self.__read_only__:
            if ignore_warning is None:
                ignore_warning = self.ignore_warnings
            if not ignore_warning:
                import warnings
                warnings.warn('Database is read only, using the cursor can '
                    + 'alter the stored data. Set self.ignore_warnings to '
                    + 'True in order to mute future warnings.', RuntimeWarning)
        return self.__connection__.cursor()

    def get_connection(self, ignore_warning=None):
        """
        Return a pysqlite connection to the database.

        You most likely want ``get_cursor()`` instead, which is used for
        executing sqlite commands on the database.

        Recommended for more advanced users only.

        EXAMPLES::

            sage: D = SQLDatabase(read_only=True)
            sage: con = D.get_connection()
            doctest:...: RuntimeWarning: Database is read only, using the connection can alter the stored data. Set self.ignore_warnings to True in order to mute future warnings.
            sage: con = D.get_connection(True)
            sage: D.ignore_warnings = True
            sage: con = D.get_connection()
            sage: t = con.execute('CREATE TABLE simon(n INTEGER, n2 INTEGER)')
            sage: for n in range(10):
            ....:   t = con.execute('INSERT INTO simon VALUES(%d,%d)'%(n,n^2))
            sage: D.show('simon')
            n                    n2
            ----------------------------------------
            0                    0
            1                    1
            2                    4
            3                    9
            4                    16
            5                    25
            6                    36
            7                    49
            8                    64
            9                    81
        """
        if self.__read_only__:
            if ignore_warning is None:
                ignore_warning = self.ignore_warnings
            if not ignore_warning:
                import warnings
                warnings.warn('Database is read only, using the connection '
                    'can alter the stored data. Set self.ignore_warnings '
                    'to True in order to mute future warnings.',
                    RuntimeWarning)
        return self.__connection__

    def create_table(self, table_name, table_skeleton):
        """
        Create a new table in the database.

        To create a table, a column structure must be specified. The form for
        this is a Python dict, for example::

            {'col1': {'sql':'INTEGER', 'index':False, 'unique':True, 'primary_key':False}, ...}

        INPUT:

        - ``table_name`` -- string
        - ``table_skeleton`` -- a double-indexed dictionary

            - outer key -- column name

                - inner key -- one of the following:

                    - ``primary_key`` -- boolean, whether column has been set
                      asprimary key
                    - ``index`` -- boolean, whether column has been set as
                      index
                    - ``unique`` -- boolean, whether column has been set as
                      unique
                    - ``sql`` -- one of ``'TEXT'``, ``'BOOLEAN'``,
                      ``'INTEGER'``, ``'REAL'``, or other user defined type

        NOTE:

        Some SQL features, such as automatically incrementing primary key,
        require the full word ``'INTEGER'``, not just ``'INT'``.

        EXAMPLES::

            sage: from collections import OrderedDict
            sage: D = SQLDatabase()
            sage: table_skeleton = OrderedDict([
            ....: ('graph6', {'sql':'TEXT', 'index':True, 'primary_key':True}),
            ....: ('vertices', {'sql':'INTEGER'}),
            ....: ('edges', {'sql':'INTEGER'})
            ....: ])
            sage: D.create_table('simon', table_skeleton)
            sage: D.show('simon')
            graph6               vertices             edges
            ------------------------------------------------------------
        """
        if self.__read_only__:
            raise RuntimeError('Cannot add table to a read only database.')
        if table_name in self.__skeleton__:
            raise ValueError('Database already has a table named %s'
                % table_name)
        if table_name.find(' ') != -1:
            raise ValueError('Table names cannot contain spaces.')
        if table_name.upper() in sqlite_keywords:
            raise ValueError('Table names cannot be a SQLite keyword.')
        create_statement = 'CREATE TABLE ' + table_name + '('
        statement = []
        index_statement = ''
        for col in table_skeleton:
            if col.find('sqlite') != -1:
                raise ValueError("Column names cannot contain 'sqlite'.")
            if col.upper() in sqlite_keywords:
                raise ValueError('Column names cannot be a SQLite keyword.')
            table_skeleton[col] = verify_column(table_skeleton[col])
            typ = table_skeleton[col]['sql']
            if verify_type(typ):
                if typ.upper() == 'NOTYPE':
                    typ = ''
                if table_skeleton[col]['primary_key']:
                    statement.append(col + ' ' + typ + ' PRIMARY KEY')
                elif table_skeleton[col]['unique']:
                    statement.append(col + ' ' + typ + ' UNIQUE')
                else:
                    statement.append(col + ' ' + typ)
                    if table_skeleton[col]['index']:
                        index_statement += 'CREATE INDEX i_%s_%s' % (table_name,
                            col) + ' ON %s(%s);\n' % (table_name, col)
        create_statement += ', '.join(statement) + ') '

        self.__connection__.execute(create_statement)
        if index_statement:
            self.__connection__.executescript(index_statement)

        self.__skeleton__[table_name] = table_skeleton

    def add_column(self, table_name, col_name, col_dict, default='NULL'):
        """
        Add a column named ``col_name`` to table ``table_name``, whose data
        types are described by ``col_dict``. The format for this is::

            {'col1':{'primary_key':False, 'unique': True, 'index':True, 'sql':'REAL'}}

        INPUT:

        - ``col_dict`` -- dictionary:

            - key -- column name

                - inner key -- one of the following:

                    - ``primary_key`` -- boolean, whether column has been set
                      as primary key
                    - ``index`` -- boolean, whether column has been set as
                      index
                    - ``unique`` -- boolean, whether column has been set as
                      unique
                    - ``sql`` -- one of ``'TEXT'``, ``'BOOLEAN'``,
                      ``'INTEGER'``, ``'REAL'``, or other user defined type

        EXAMPLES::

            sage: from collections import OrderedDict
            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', OrderedDict([('n', {'sql':'INTEGER', 'index':True})]))
            sage: for n in range(20): MonicPolys.add_row('simon', (n,))
            sage: MonicPolys.add_column('simon', 'n_squared', {'sql':'INTEGER', 'index':False}, 0)
            sage: MonicPolys.show('simon')
            n                    n_squared
            ----------------------------------------
            0                    0
            1                    0
            2                    0
            3                    0
            4                    0
            5                    0
            6                    0
            7                    0
            8                    0
            9                    0
            10                   0
            11                   0
            12                   0
            13                   0
            14                   0
            15                   0
            16                   0
            17                   0
            18                   0
            19                   0
        """
        if self.__read_only__:
            raise RuntimeError('Cannot add columns to a read only database.')
        # Check input:
        if col_name.find('sqlite') != -1:
            raise ValueError("Column names cannot contain 'sqlite'.")
        if col_name.upper() in sqlite_keywords:
            raise ValueError("Column names cannot be SQLite keywords.")
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        if col_name in self.__skeleton__[table_name]:
            raise ValueError("Table %s already has column %s." % (table_name,col_name))

        # Update the skeleton:
        self.__skeleton__[table_name][col_name] = verify_column(col_dict)

        try:
            self._rebuild_table(table_name, col_name, default)
        except sqlite.Error as e:
            # delete added column from skeleton
            self.__skeleton__[table_name].pop(col_name)

            print('A sqlite error occurred: ', e.args[0])

    def _rebuild_table(self, table_name, col_name=None, default=''):
        """
        Rebuild the table ``table_name`` adding column ``col_name`` if not
        ``None``. If a new column is added, each rows' value is set to
        ``default``.

        Used in the methods ``add_column``, ``drop_column``, ``make_unique``,
        ``drop_unique``, ``make_primary_key``, and ``drop_primary_key`` to get
        around the limitations of SQLite's ``ALTER TABLE``. We have set up the
        creation of a temporary database in this method. Please feel free to
        improve on this by sending a patch or suggestion.

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}})
            sage: for n in range(20): MonicPolys.add_row('simon', (n,))
            sage: MonicPolys.show('simon')
            n
            --------------------
            0
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
            13
            14
            15
            16
            17
            18
            19
            sage: MonicPolys._rebuild_table('simon')
            sage: MonicPolys.show('simon')
            n
            --------------------
            0
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
            13
            14
            15
            16
            17
            18
            19
        """
        cols = []
        cols_attr = []
        table_skeleton = self.__skeleton__[table_name]
        # gather column names and attributes for new table
        for col in table_skeleton:
            cols.append(col)
            attr_str = col + ' ' + table_skeleton[col]['sql']
            if table_skeleton[col]['primary_key']:
                attr_str += ' PRIMARY KEY'
            elif table_skeleton[col]['unique']:
                attr_str += ' UNIQUE'
            cols_attr.append(attr_str)

        original = list(cols)

        if col_name is not None:
            original[original.index(col_name)] = str(default)

        original = ', '.join(original)
        cols = ', '.join(cols)
        cols_attr = ', '.join(cols_attr)

        # Silly SQLite -- we have to make a temp table to hold info...
        self.__connection__.executescript("""
            CREATE TEMPORARY TABLE spam(%s);
            INSERT INTO spam SELECT %s FROM %s;
            DROP TABLE %s;
            CREATE TABLE %s (%s);
            """ % (cols_attr, original, table_name, table_name, table_name, cols_attr))
        # Update indices in new table
        skeleton = self.__skeleton__[table_name]
        index_statement = ''.join(f'CREATE INDEX i_{table_name}_{col} ON '
                                  + f'{table_name}({col});\n'
                                  for col in skeleton
                                  if skeleton[col]['index']
                                  and not skeleton[col]['primary_key'])
        if index_statement:
            self.__connection__.executescript(index_statement)

        # Now we can plop our data into the *new* table:
        self.__connection__.executescript("""
            INSERT INTO %s SELECT %s FROM spam;
            DROP TABLE spam;
            """ % (table_name, cols))

        self.vacuum()

    def drop_column(self, table_name, col_name):
        """
        Drop the column ``col_name`` from table ``table_name``.

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}})
            sage: for n in range(20): MonicPolys.add_row('simon', (n,))
            sage: MonicPolys.add_column('simon', 'n_squared', {'sql':'INTEGER'}, 0)
            sage: MonicPolys.drop_column('simon', 'n_squared')
            sage: MonicPolys.show('simon')
            n
            --------------------
            0
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
            13
            14
            15
            16
            17
            18
            19
        """
        if self.__read_only__:
            raise RuntimeError('Cannot drop columns in a read only database.')
        # Check input:
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        if col_name not in self.__skeleton__[table_name]:
            raise ValueError("Table %s has no column %s." % (table_name,col_name))

        # Update the skeleton:
        self.__skeleton__[table_name].pop(col_name)

        self._rebuild_table(table_name)

    def rename_table(self, table_name, new_name):
        """
        Rename the table ``table_name`` to ``new_name``.

        EXAMPLES::

            sage: D = SQLDatabase()
            sage: D.create_table('simon',{'col1':{'sql':'INTEGER'}})
            sage: D.show('simon')
            col1
            --------------------
            sage: D.rename_table('simon', 'lucy')
            sage: D.show('simon')
            Traceback (most recent call last):
            ...
            RuntimeError: Failure to fetch data.
            sage: D.show('lucy')
            col1
            --------------------
        """
        if self.__read_only__:
            raise RuntimeError('Cannot rename tables in a read only database.')
        # Check input:
        if table_name not in self.__skeleton__:
            raise ValueError('Database has no table %s.' % table_name)
        if new_name in self.__skeleton__:
            raise ValueError('Database already has table %s.' % new_name)

        self.__connection__.execute('ALTER TABLE %s RENAME TO ' % table_name
            + new_name)

        # Update skeleton:
        self.__skeleton__[new_name] = self.__skeleton__.pop(table_name)

    def drop_table(self, table_name):
        """
        Delete table ``table_name`` from database.

        INPUT:

        - ``table_name`` -- string

        EXAMPLES::

            sage: D = SQLDatabase()
            sage: D.create_table('simon',{'col1':{'sql':'INTEGER'}})
            sage: D.show('simon')
            col1
            --------------------
            sage: D.drop_table('simon')
            sage: D.get_skeleton()
            {}
        """
        if self.__read_only__:
            raise RuntimeError('Cannot drop tables from a read only '
                + 'database.')
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)

        self.__connection__.execute('DROP TABLE ' + table_name)

        # Update Skeleton
        self.__skeleton__.pop(table_name)

    def drop_data_from_table(self, table_name):
        """
        Remove all rows from ``table_name``.

        EXAMPLES::

            sage: D = SQLDatabase()
            sage: D.create_table('simon',{'col1':{'sql':'INTEGER'}})
            sage: D.add_row('simon',(9,))
            sage: D.show('simon')
            col1
            --------------------
            9
            sage: D.drop_data_from_table('simon')
            sage: D.show('simon')
            col1
            --------------------
        """
        if self.__read_only__:
            raise RuntimeError('Cannot remove data from a read only database.')
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        self.__connection__.execute('DELETE FROM ' + table_name)

    def make_index(self, col_name, table_name, unique=False):
        """
        Set the column ``col_name`` in table ``table_name`` to be an index,
        that is, a column set up to do quick searches on.

        INPUT:

        - ``col_name`` -- string
        - ``table_name`` -- string
        - ``unique`` -- requires that there are no multiple entries in the
          column, makes searching faster

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: MonicPolys.make_index('n2','simon')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'n2': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False}}}
        """
        if self.__read_only__:
            raise RuntimeError('Cannot modify a read only database.')
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        if col_name not in self.__skeleton__[table_name]:
            raise ValueError("Table %s has no column %s." % (table_name,col_name))

        if unique:
            index_string = 'CREATE UNIQUE INDEX ' + col_name + ' ON ' \
                + table_name + ' (' + col_name + ')'
        else:
            index_string = 'CREATE INDEX ' + col_name + ' ON ' + table_name \
                + ' (' + col_name + ')'
        cur = self.__connection__.cursor()
        cur.execute(index_string)

        # Update Skeleton
        self.__skeleton__[table_name][col_name]['index'] = True
        if unique:
            self.__skeleton[table_name][col_name]['unique'] = True

    def drop_index(self, table_name, index_name):
        """
        Set the column ``index_name`` in table ``table_name`` to not be an
        index. See ``make_index()``

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: MonicPolys.drop_index('simon', 'n')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'n2': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False}}}
        """
        if self.__read_only__:
            raise RuntimeError('Cannot modify a read only database.')
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        if index_name not in self.__skeleton__[table_name]:
            raise ValueError("Table %s has no column %s." % (table_name,
                                                             index_name))
        if not self.__skeleton__[table_name][index_name]['index']:
            return # silently

        cur = self.__connection__.cursor()
        cur.execute('DROP INDEX i_' + table_name + '_' + index_name)

        # Update Skeleton
        self.__skeleton__[table_name][index_name]['index'] = False

    def make_unique(self, table_name, col_name):
        """
        Set the column ``col_name`` in table ``table_name`` to store unique
        values.

        NOTE:

        This function only adds the requirement for entries in ``col_name`` to
        be unique, it does not change the values.

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: MonicPolys.make_unique('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'n2': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': True}}}
        """
        if self.__read_only__:
            raise RuntimeError('Cannot modify a read only database')
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        if col_name not in self.__skeleton__[table_name]:
            raise ValueError("Table %s has no column %s." % (table_name, col_name))

        # Update the skeleton:
        self.__skeleton__[table_name][col_name]['unique'] = True

        self._rebuild_table(table_name)

    def drop_unique(self, table_name, col_name):
        """
        Set the column ``col_name`` in table ``table_name`` not store unique
        values.

        NOTE:

        This function only removes the requirement for entries in ``col_name``
        to be unique, it does not delete it.

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: MonicPolys.make_unique('simon', 'n2')
            sage: MonicPolys.drop_unique('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'n2': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False}}}
        """
        if self.__read_only__:
            raise RuntimeError('Cannot modify a read only database.')
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        if col_name not in self.__skeleton__[table_name]:
            raise ValueError("Table %s has no column %s." % (table_name, col_name))
        if self.__skeleton__[table_name][col_name]['primary_key']:
            raise ValueError("Primary keys must be unique.")

        # Update the skeleton:
        self.__skeleton__[table_name][col_name]['unique'] = False

        self._rebuild_table(table_name)

    def make_primary_key(self, table_name, col_name):
        """
        Set the column ``col_name`` in table ``table_name`` to be a primary key.

        A primary key is something like an index, but its main purpose is to
        link different tables together. This allows searches to be executed on
        multiple tables that represent maybe different data about the same
        objects.

        NOTE:

        Some SQL features, such as automatically incrementing primary key,
        require the full word ``'INTEGER'``, not just ``'INT'``.

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: MonicPolys.make_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'n2': {'index': False,
               'primary_key': True,
               'sql': 'INTEGER',
               'unique': True}}}
        """
        if self.__read_only__:
            raise RuntimeError('Cannot modify a read only database.')
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        if col_name not in self.__skeleton__[table_name]:
            raise ValueError("Table %s has no column %s." % (table_name, col_name))

        # Update the skeleton:
        self.__skeleton__[table_name][col_name]['primary_key'] = True
        self.__skeleton__[table_name][col_name]['unique'] = True

        self._rebuild_table(table_name)

    def drop_primary_key(self, table_name, col_name):
        """
        Set the column ``col_name`` in table ``table_name`` not to be a primary
        key.

        A primary key is something like an index, but its main purpose is to
        link different tables together. This allows searches to be executed on
        multiple tables that represent maybe different data about the same
        objects.

        NOTE:

        This function only changes the column to be non-primary, it does not
        delete it.

        EXAMPLES::

            sage: MonicPolys = SQLDatabase()
            sage: MonicPolys.create_table('simon', {'n':{'sql':'INTEGER', 'index':True}, 'n2':{'sql':'INTEGER'}})
            sage: MonicPolys.make_primary_key('simon', 'n2')
            sage: MonicPolys.drop_primary_key('simon', 'n2')
            sage: MonicPolys.get_skeleton()
            {'simon': {'n': {'index': True,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': False},
              'n2': {'index': False,
               'primary_key': False,
               'sql': 'INTEGER',
               'unique': True}}}
        """
        if self.__read_only__:
            raise RuntimeError('Cannot modify a read only database.')
        if table_name not in self.__skeleton__:
            raise ValueError("Database has no table %s." % table_name)
        if col_name not in self.__skeleton__[table_name]:
            raise ValueError("Table %s has no column %s." % (table_name,col_name))
        if not self.__skeleton__[table_name][col_name]['primary_key']:
            return # silently

        # Update the skeleton:
        self.__skeleton__[table_name][col_name]['primary_key'] = False

        self._rebuild_table(table_name)

    def add_row(self, table_name, values, entry_order=None):
        """
        Add the row described by ``values`` to the table ``table_name``. Values
        should be a tuple, of same length and order as columns in given table.

        NOTE:

        If ``values`` is of length one, be sure to specify that it is a tuple of
        length one, by using a comma, e.g.::

            sage: values = (6,)

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_row('simon',(0,1))
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(0, 1)]
        """
        self.add_rows(table_name, [values], entry_order)

    def delete_rows(self, query):
        """
        Use a ``SQLQuery`` instance to modify (delete rows from) the
        database.

        ``SQLQuery`` must have no join statements.  (As of now, you can only
        delete from one table at a time -- ideas and patches are welcome).

        To remove all data that satisfies a ``SQLQuery``, send the query as an
        argument to ``delete_rows``.  Be careful, test your query first.

        Recommended use:  have some kind of row identification primary
        key column that you use as a parameter in the query.  (See example
        below).

        INPUT:

        - ``query`` -- a ``SQLQuery`` (Delete the rows returned when query is
          run)

        EXAMPLES::

            sage: from collections import OrderedDict
            sage: DB = SQLDatabase()
            sage: DB.create_table('lucy', OrderedDict([
            ....: ('id', {'sql':'INTEGER', 'primary_key':True, 'index':True}),
            ....: ('a1', {'sql':'bool'}),
            ....: ('b2', {'sql':'int'})]))
            sage: DB.add_rows('lucy', [(0,1,1),(1,1,4),(2,0,7),(3,1,384), (4,1,978932)],['id','a1','b2'])
            sage: DB.show('lucy')
            id                   a1                   b2
            ------------------------------------------------------------
            0                    1                    1
            1                    1                    4
            2                    0                    7
            3                    1                    384
            4                    1                    978932
            sage: Q = SQLQuery(DB, {'table_name':'lucy', 'display_cols':['id','a1','b2'], 'expression':['id','>=',3]})
            sage: DB.delete_rows(Q)
            sage: DB.show('lucy')
            id                   a1                   b2
            ------------------------------------------------------------
            0                    1                    1
            1                    1                    4
            2                    0                    7
        """
        if self.__read_only__:
            raise RuntimeError('Cannot delete rows from a read only database.')
        # Check query is associated with this database
        if not isinstance(query, SQLQuery):
            raise TypeError('%s is not a valid SQLQuery' % query)
        if query.__database__ is not self:
            raise ValueError('%s is not associated to this database.' % query)
        if (query.__query_string__).find(' JOIN ') != -1:
            raise ValueError(f'{query} is not a valid query. Can only '
                             'delete from one table at a time.')

        delete_statement = re.sub('SELECT .* FROM', 'DELETE FROM',
                                  query.__query_string__)

        try:
            cur = self.get_cursor()
            cur.execute(delete_statement, query.__param_tuple__)
        except Exception:
            raise RuntimeError('Failure to complete delete. Check your data.')

    def add_rows(self, table_name, rows, entry_order=None):
        """
        INPUT:

        - ``rows`` -- list of tuples that represent one row of data to add
          (types should match col types in order)
        - ``entry_order`` -- an ordered list or tuple overrides normal order
          with user defined order

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_rows('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_rows('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: cur = DB.get_cursor()
            sage: (cur.execute('select * from simon')).fetchall()
            [(0, 0), (1, 1), (1, 2), (0, 0), (0, 4), (1, 5)]
        """
        if self.__read_only__:
            raise RuntimeError('Cannot add rows to read only database.')
        quest = '(' + ', '.join('?' for i in rows[0]) + ')'
        strows = [tuple(str(entry) for entry in row) for row in rows]

        if entry_order is not None:
            self.__connection__.executemany('INSERT INTO ' + table_name
                + str(tuple(entry_order)) + ' VALUES ' + quest, strows)
        else:
            self.__connection__.executemany('INSERT INTO ' + table_name
                + ' VALUES ' + quest, strows)

    add_data = add_rows

    def vacuum(self):
        """
        Clean the extra hard disk space used up by a database that has
        recently shrunk.

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_row('simon',(0,1))
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_data('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: DB.drop_column('simon','b2')
            sage: DB.commit()
            sage: DB.vacuum()
        """
        self.__connection__.execute('VACUUM')

    def commit(self):
        """
        Commit changes to file.

        EXAMPLES::

            sage: DB = SQLDatabase()
            sage: DB.create_table('simon',{'a1':{'sql':'bool'}, 'b2':{'sql':'int'}})
            sage: DB.add_row('simon',(0,1))
            sage: DB.add_data('simon',[(0,0),(1,1),(1,2)])
            sage: DB.add_data('simon',[(0,0),(4,0),(5,1)], ['b2','a1'])
            sage: DB.drop_column('simon','b2')
            sage: DB.commit()
            sage: DB.vacuum()
        """
        if self.__read_only__:
            raise RuntimeError("Cannot commit read only database.")
        try:
            self.__connection__.execute('COMMIT')
        except sqlite.OperationalError:
            # Not sure why this throws an exception - but without it,
            #       the changes are not committed so it is necessary.
            pass
