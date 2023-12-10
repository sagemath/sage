# sage_setup: distribution = sagemath-graphs
r"""
Lovász theta-function of graphs

AUTHORS:

- Dima Pasechnik (2015-06-30): Initial version

REFERENCE:

[Lov1979]_

Functions
---------
"""


def lovasz_theta(graph):
    r"""
    Return the value of Lovász theta-function of graph.

    For a graph `G` this function is denoted by `\theta(G)`, and it can be
    computed in polynomial time. Mathematically, its most important property is
    the following:

    .. MATH::

        \alpha(G)\leq\theta(G)\leq\chi(\overline{G})

    with `\alpha(G)` and `\chi(\overline{G})` being, respectively, the maximum
    size of an :meth:`independent set <sage.graphs.graph.Graph.independent_set>`
    set of `G` and the :meth:`chromatic number
    <sage.graphs.graph.Graph.chromatic_number>` of the :meth:`complement
    <sage.graphs.generic_graph.GenericGraph.complement>` `\overline{G}` of `G`.

    For more information, see the :wikipedia:`Lovász_number`.

    .. NOTE::

        - Implemented for undirected graphs only. Use ``to_undirected``
          to convert a digraph to an undirected graph.

        - This function requires the optional package ``csdp``, which you can
          install with ``sage -i csdp``.

    EXAMPLES::

          sage: C = graphs.PetersenGraph()
          sage: C.lovasz_theta()                             # optional - csdp
          4.0
          sage: graphs.CycleGraph(5).lovasz_theta()          # optional - csdp
          2.236068

    TESTS::

        sage: g = Graph()
        sage: g.lovasz_theta() # indirect doctest
        0
    """
    n = graph.order()
    if not n:
        return 0

    from networkx import write_edgelist
    from sage.misc.temporary_file import tmp_filename
    import subprocess

    from sage.features.csdp import CSDP
    CSDP().require()

    g = graph.relabel(inplace=False, perm=range(1, n + 1)).networkx_graph()
    tf_name = tmp_filename()
    with open(tf_name, 'wb') as tf:
        tf.write("{}\n{}\n".format(n, g.number_of_edges()).encode())
        write_edgelist(g, tf, data=False)
    lines = subprocess.check_output(['theta', tf_name])
    return float(lines.split()[-1])
