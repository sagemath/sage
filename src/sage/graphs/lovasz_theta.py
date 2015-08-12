# -*- coding: utf-8 -*-
r"""
Lovász theta-function of graphs

AUTHORS:

- Dima Pasechnik (2015-06-30): Initial version

REFERENCE:

.. [Lovasz1979] László Lovász,
  "On the Shannon capacity of a graph",
  IEEE Trans. Inf. Th. 25(1979), 1-7.

Functions
---------
"""

def lovasz_theta(graph):
    r"""
    Return the value of Lovász theta-function of graph

    For a graph `G` this function is denoted by `\theta(G)`, and it can be
    computed in polynomial time. Mathematically, its most important property is the following:

    .. MATH::

        \alpha(G)\leq\theta(G)\leq\chi(\overline{G})

    with `\alpha(G)` and `\chi(\overline{G})` being, respectively, the maximum
    size of an :meth:`independent set <sage.graphs.graph.Graph.independent_set>`
    set of `G` and the :meth:`chromatic number
    <sage.graphs.graph.Graph.chromatic_number>` of the :meth:`complement
    <sage.graphs.generic_graph.GenericGraph.complement>` `\overline{G}` of `G`.

    For more information, see the :wikipedia:`Lovász_number`.

    .. NOTE::

        - Implemented for undirected graphs only. Use to_undirected to convert a
          digraph to an undirected graph.

        - This function requires the optional package ``csdp``, which you can
          install with with ``sage -i csdp``.

    EXAMPLES::

          sage: C=graphs.PetersenGraph()
          sage: C.lovasz_theta()                             # optional csdp
          4.0
          sage: graphs.CycleGraph(5).lovasz_theta()          # optional csdp
          2.236068

    TEST::

        sage: g = Graph()
        sage: g.lovasz_theta() # indirect doctest
        0
    """
    n = graph.order()
    if n == 0:
        return 0

    from networkx import write_edgelist
    from sage.misc.temporary_file import tmp_filename
    import os, subprocess
    from sage.env import SAGE_LOCAL
    from sage.misc.package import is_package_installed

    if not is_package_installed('csdp'):
        raise NotImplementedError("Package csdp is required. Please install it with 'sage -i csdp'.")

    g = graph.relabel(inplace=False, perm=range(1,n+1)).networkx_graph()
    tf_name = tmp_filename()
    tf = open(tf_name, 'wb')
    tf.write(str(n)+'\n'+str(g.number_of_edges())+'\n')
    write_edgelist(g, tf, data=False)
    tf.close()
    lines = subprocess.check_output([os.path.join(SAGE_LOCAL, 'bin', 'theta'), tf_name])
    return float(lines.split()[-1])
