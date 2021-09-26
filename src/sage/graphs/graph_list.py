r"""
Lists of graphs

AUTHORS:

- Robert L. Miller (2007-02-10): initial version

- Emily A. Kirkman (2007-02-13): added show functions
  (to_graphics_array and show_graphs)
"""

# ****************************************************************************
#           Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.misc import try_read


def from_whatever(data):
    r"""
    Return a list of Sage Graphs, given a list of whatever kind of data.

    INPUT:

    - ``data`` -- can be a string, a list/iterable of strings, or a readable
      file-like object

    EXAMPLES::

        sage: l = ['N@@?N@UGAGG?gGlKCMO', ':P_`cBaC_ACd`C_@BC`ABDHaEH_@BF_@CHIK_@BCEHKL_BIKM_BFGHI']
        sage: graphs_list.from_whatever(l)
        [Graph on 15 vertices, Looped multi-graph on 17 vertices]
        sage: graphs_list.from_whatever('\n'.join(l))
        [Graph on 15 vertices, Looped multi-graph on 17 vertices]

    This example happens to be a mix a sparse and non-sparse graphs, so we don't
    explicitly put a ``.g6`` or ``.s6`` extension, which implies just one or the
    other::

        sage: filename = tmp_filename()
        sage: with open(filename, 'w') as fobj:
        ....:     _ = fobj.write('\n'.join(l))
        sage: with open(filename) as fobj:
        ....:     graphs_list.from_whatever(fobj)
        [Graph on 15 vertices, Looped multi-graph on 17 vertices]
    """
    return _from_whatever(data)


def _from_whatever(data, fmt=None):
    """
    Implementation details of :func:`from_whatever`.

    INPUT:

    - ``data`` -- can be a string, a list/iterable of strings, or a readable
      file-like object

    - ``fmt`` -- string (default: ``None``); format of ``data``. It can be
      either ``'graph6'``, ``sparse6``, or ``None``, with the latter case
      indicating that the ``Graph`` constructor should determine this for
      itself

    EXAMPLES::

        sage: l = ['N@@?N@UGAGG?gGlKCMO', ':P_`cBaC_ACd`C_@BC`ABDHaEH_@BF_@CHIK_@BCEHKL_BIKM_BFGHI']
        sage: graphs_list.from_whatever(l)
        [Graph on 15 vertices, Looped multi-graph on 17 vertices]
    """
    from sage.graphs.graph import Graph

    if isinstance(data, str):
        lines = data.splitlines()
    else:
        lines = try_read(data, splitlines=True)

        if lines is not None and fmt is None:
            # In this case the format should be 'forced' by the filename
            if hasattr(data, 'name'):
                if data.name.endswith('.g6'):
                    fmt = 'graph6'
                elif data.name.endswith('.s6'):
                    fmt = 'sparse6'
        else:
            try:
                lines = iter(data)
            except TypeError:
                raise TypeError(
                    "must be a string, an iterable of strings, or a readable "
                    "file-like object")

    if fmt == 'graph6':
        kwargs = {'format': fmt}
    elif fmt == 'sparse6':
        kwargs = {'format': fmt, 'sparse': True}  # probably implied?
    else:
        kwargs = {}  # We let Graph guess

    out = []
    for line in lines:
        if not isinstance(line, str):
            raise TypeError("must be an iterable of strings")
        line = line.strip()
        if not line:
            continue

        if '\n' in line:
            out.append(_from_whatever(line.splitlines(), fmt=fmt))
        else:
            out.append(Graph(line, **kwargs))

    return out


def from_graph6(data):
    """
    Return a list of Sage Graphs, given a list of graph6 data.

    INPUT:

    - ``data`` -- can be a string, a list of strings, or a file stream

    EXAMPLES::

        sage: l = ['N@@?N@UGAGG?gGlKCMO', 'XsGGWOW?CC?C@HQKHqOjYKC_uHWGX?P?~TqIKA`OA@SAOEcEA??']
        sage: graphs_list.from_graph6(l)
        [Graph on 15 vertices, Graph on 25 vertices]
    """
    return _from_whatever(data, fmt='graph6')


def from_sparse6(data):
    """
    Return a list of Sage Graphs, given a list of sparse6 data.

    INPUT:

    - ``data`` -- can be a string, a list of strings, or a file stream

    EXAMPLES::

        sage: l = [':P_`cBaC_ACd`C_@BC`ABDHaEH_@BF_@CHIK_@BCEHKL_BIKM_BFGHI', ':f`??KO?B_OOSCGE_?OWONDBO?GOJBDO?_SSJdApcOIG`?og_UKEbg?_SKFq@[CCBA`p?oYMFp@gw]Qaa@xEMHDb@hMCBCbQ@ECHEcAKKQKFPOwo[PIDQ{KIHEcQPOkVKEW_WMNKqPWwcRKOOWSKIGCqhWt??___WMJFCahWzEBa`xOu[MpPPKqYNoOOOKHHDBPs|??__gWMKEcAHKgTLErqA?A@a@G{kVLErs?GDBA@XCs\\NggWSOJIDbHh@?A@aF']
        sage: graphs_list.from_sparse6(l)
        [Looped multi-graph on 17 vertices, Looped multi-graph on 39 vertices]
    """
    return _from_whatever(data, fmt='sparse6')


def to_graph6(graphs, file=None, output_list=False):
    r"""
    Convert a list of Sage graphs to a single string of graph6 graphs.

    If ``file`` is specified, then the string will be written quietly to the
    file.  If ``output_list`` is ``True``, then a list of strings will be
    returned, one string per graph.

    INPUT:

    - ``graphs`` -- a Python list of Sage Graphs

    - ``file`` -- (optional) a file stream to write to (must be in 'w' mode)

    - ``output_list`` -- boolean (default: ``False``); whether to return a
      string (when set to ``True``) or a list of strings. This parameter is
      ignored if file gets specified

    EXAMPLES::

        sage: l = [graphs.DodecahedralGraph(), graphs.PetersenGraph()]
        sage: graphs_list.to_graph6(l)
        'ShCHGD@?K?_@?@?C_GGG@??cG?G?GK_?C\nIheA@GUAo\n'
    """
    return _to_graph6(graphs, file=file, output_list=output_list)


def to_sparse6(graphs, file=None, output_list=False):
    r"""
    Convert a list of Sage graphs to a single string of sparse6 graphs.

    If ``file`` is specified, then the string will be written quietly to the
    file.  If ``output_list`` is ``True``, then a list of strings will be
    returned, one string per graph.

    INPUT:

    - ``graphs`` -- a Python list of Sage Graphs

    - ``file`` -- (optional) a file stream to write to (must be in 'w' mode)

    - ``output_list`` -- boolean (default: ``False``); whether to return a
      string (when set to ``True``) or a list of strings. This parameter is
      ignored if file gets specified

    EXAMPLES::

        sage: l = [graphs.DodecahedralGraph(), graphs.PetersenGraph()]
        sage: graphs_list.to_sparse6(l)
        ':S_`abcaDe`Fg_HijhKfLdMkNcOjP_BQ\n:I`ES@obGkqegW~\n'
    """
    return _to_graph6(graphs, file=file, output_list=output_list, sparse=True)


def _to_graph6(graphs, file=None, output_list=False, sparse=False):
    r"""
    Internal implementation of :func:`to_graph6` and :func:`to_sparse6`.

    EXAMPLES::

        sage: l = [graphs.DodecahedralGraph(), graphs.PetersenGraph()]
        sage: graphs_list._to_graph6(l, sparse=False)
        'ShCHGD@?K?_@?@?C_GGG@??cG?G?GK_?C\nIheA@GUAo\n'
        sage: graphs_list._to_graph6(l, sparse=True)
        ':S_`abcaDe`Fg_HijhKfLdMkNcOjP_BQ\n:I`ES@obGkqegW~\n'
    """
    if sparse:
        method = 'sparse6_string'
    else:
        method = 'graph6_string'

    strs = [getattr(g, method)() for g in graphs]

    if file or not output_list:
        strs = '\n'.join(strs) + '\n'

    if file is None:
        return strs

    file.write(strs)
    file.flush()


def to_graphics_array(graph_list, **kwds):
    """
    Draw all graphs in a graphics array

    INPUT:

    - ``graph_list`` -- a Python list of Sage Graphs

    GRAPH PLOTTING:

    Defaults to circular layout for graphs. This allows for a nicer display in a
    small area and takes much less time to compute than the spring- layout
    algorithm for many graphs.

    EXAMPLES::

        sage: glist = []
        sage: for i in range(999):
        ....:     glist.append(graphs.RandomGNP(6, .45))
        sage: garray = graphs_list.to_graphics_array(glist)
        sage: garray.nrows(), garray.ncols()
        (250, 4)

    See the .plot() or .show() documentation for an individual graph for
    options, all of which are available from :func:`to_graphics_array`::

        sage: glist = []
        sage: for _ in range(10):
        ....:     glist.append(graphs.RandomLobster(41, .3, .4))
        sage: graphs_list.to_graphics_array(glist, layout='spring', vertex_size=20)
        Graphics Array of size 3 x 4
    """
    from sage.graphs import graph
    plist = []
    for graph_i in graph_list:
        if isinstance(graph_i, graph.GenericGraph):
            pos = graph_i.get_pos()
            if pos is None:
                if 'layout' not in kwds:
                    kwds['layout'] = 'circular'
                if 'vertex_size' not in kwds:
                    kwds['vertex_size'] = 50
                if 'vertex_labels' not in kwds:
                    kwds['vertex_labels'] = False
                kwds['graph_border'] = True
                plist.append(graph_i.plot(**kwds))
            else:
                plist.append(graph_i.plot(pos=pos, vertex_size=50,
                                          vertex_labels=False,
                                          graph_border=True))
        else:
            raise TypeError('param list must be a list of Sage (di)graphs.')
    from sage.plot.plot import graphics_array
    return graphics_array(plist, ncols=4)


def show_graphs(graph_list, **kwds):
    """
    Show a maximum of 20 graphs from ``graph_list`` in a sage graphics array.

    If more than 20 graphs are given in the list argument, then it will display
    one graphics array after another with each containing at most 20 graphs.

    Note that to save the image output from the notebook, you must save each
    graphics array individually. (There will be a small space between graphics
    arrays).

    INPUT:

    - ``graph_list`` -- a Python list of Sage Graphs

    GRAPH PLOTTING: Defaults to circular layout for graphs. This allows for a
    nicer display in a small area and takes much less time to compute than the
    spring-layout algorithm for many graphs.

    EXAMPLES: Create a list of graphs::

        sage: glist = []
        sage: glist.append(graphs.CompleteGraph(6))
        sage: glist.append(graphs.CompleteBipartiteGraph(4, 5))
        sage: glist.append(graphs.BarbellGraph(7, 4))
        sage: glist.append(graphs.CycleGraph(15))
        sage: glist.append(graphs.DiamondGraph())
        sage: glist.append(graphs.GemGraph())
        sage: glist.append(graphs.DartGraph())
        sage: glist.append(graphs.ForkGraph())
        sage: glist.append(graphs.HouseGraph())
        sage: glist.append(graphs.HouseXGraph())
        sage: glist.append(graphs.KrackhardtKiteGraph())
        sage: glist.append(graphs.LadderGraph(5))
        sage: glist.append(graphs.LollipopGraph(5, 6))
        sage: glist.append(graphs.PathGraph(15))
        sage: glist.append(graphs.PetersenGraph())
        sage: glist.append(graphs.StarGraph(17))
        sage: glist.append(graphs.WheelGraph(9))

    Check that length is <= 20::

        sage: len(glist)
        17

    Show the graphs in a graphics array::

        sage: graphs_list.show_graphs(glist)

    Example where more than one graphics array is used::

        sage: gq = GraphQuery(display_cols=['graph6'], num_vertices=5)
        sage: g = gq.get_graphs_list()
        sage: len(g)
        34
        sage: graphs_list.show_graphs(g)

    See the .plot() or .show() documentation for an individual graph for
    options, all of which are available from :func:`to_graphics_array`::

        sage: glist = []
        sage: for _ in range(10):
        ....:     glist.append(graphs.RandomLobster(41, .3, .4))
        sage: graphs_list.show_graphs(glist, layout='spring', vertex_size=20)
    """
    graph_list = list(graph_list)
    for i in range(len(graph_list) // 20 + 1):
        graph_slice = graph_list[20 * i: 20 * (i + 1)]
        to_graphics_array(graph_slice, **kwds).show()
