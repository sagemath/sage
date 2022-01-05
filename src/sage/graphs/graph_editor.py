r"""
Graph editor
"""
# ****************************************************************************
#      Copyright (C) 2009   Radoslav Kirov
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import sys

from .graph_generators import graphs
from sage.misc.html import html


def graph_to_js(g):
    """
    Returns a string representation of a :class:`Graph` instance
    usable by the :func:`graph_editor`.  The encoded information is
    the number of vertices, their 2D positions, and a list of edges.

    INPUT:

    - ``g`` - a :class:`Graph` instance

    OUTPUT:

    - a string

    EXAMPLES::

        sage: from sage.graphs.graph_editor import graph_to_js
        sage: G = graphs.CompleteGraph(4)
        sage: graph_to_js(G)
        'num_vertices=4;edges=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]];pos=[[0.5,0.0],[0.0,0.5],[0.5,1.0],[1.0,0.5]];'
        sage: graph_to_js(graphs.StarGraph(2))
        'num_vertices=3;edges=[[0,1],[0,2]];pos=[[0.0,0.5],[0.0,0.0],[0.0,1.0]];'
    """
    string = ''
    vertex_list = list(g.get_vertices())
    string += 'num_vertices=' + str(len(vertex_list)) + ';'
    string += 'edges=['
    for i, e in enumerate(g.edges()):
        if i:
            string += ','
        string += '[' + str(vertex_list.index(e[0])) + ',' + str(vertex_list.index(e[1])) + ']'
    string += '];'
    string += 'pos=['
    pos = g.get_pos()
    max_x = max([i[0] for i in pos.values()])
    max_y = max([i[1] for i in pos.values()])
    min_x = min([i[0] for i in pos.values()])
    min_y = min([i[1] for i in pos.values()])
    if max_x == 0:
        max_x = 1
    if max_y == 0:
        max_y = 1
    for i, v in enumerate(vertex_list):
        if i:
            string += ','
        new_pos = [float(pos[v][0] - min_x) / (max_x - min_x),
                   1.0 - float(pos[v][1] - min_y) / (max_y - min_y)]
        string += str(new_pos)
    string += '];'
    string = string.replace(' ', '')
    return string


def graph_editor(graph=None, graph_name=None,
                 replace_input=True, **layout_options):
    """
    Opens a graph editor in the Sage notebook.

    INPUT:

    - ``graph`` - a :class:`Graph` instance (default:
      graphs.CompleteGraph(2)); the graph to edit

    - ``graph_name`` - a string (default: None); the variable name to
      use for the updated instance; by default, this function attempts
      to determine the name automatically

    - ``replace_input`` - a boolean (default: True); whether to
      replace the text in the input cell with the updated graph data
      when "Save" is clicked; if this is False, the data is **still**
      evaluated as if it had been entered in the cell

    EXAMPLES::

        sage: g = graphs.CompleteGraph(3)
        sage: graph_editor(g)                       # not tested
        sage: graph_editor(graphs.HouseGraph())     # not tested
        sage: graph_editor(graph_name='my_graph')   # not tested
        sage: h = graphs.StarGraph(6)
        sage: graph_editor(h, replace_input=False)  # not tested
    """
    import sagenb.notebook.interact
    if graph is None:
        graph = graphs.CompleteGraph(2)

    return "This graph editor only runs in the deprecated Sage notebook."

    graph.layout(save_pos=True, **layout_options)

    if graph_name is None:
        graph_name = ''
        locs = sys._getframe(1).f_locals
        for var in locs:
            if id(locs[var]) == id(graph):
                graph_name = var

    cell_id = sagenb.notebook.interact.SAGE_CELL_ID

    # TODO: Put reasonable checks for large graphs, before disaster
    # occurs (i.e., breaks browser).

    close_button = r"""<button onclick="cell_delete_output(%(cell_id)s);">Close</button>""" % locals()

    if replace_input:
        eval_strategy = r"""
    f += ' graph_editor(' + g[2] + ');'
    \$('#cell_input_%(cell_id)s').val(f);
    cell_input_resize(%(cell_id)s);
    evaluate_cell(%(cell_id)s, false);
""" % locals()
    else:
        eval_strategy = r"""
    saved_input = \$('#cell_input_%(cell_id)s').val();
    \$('#cell_input_%(cell_id)s').val(f);
    evaluate_cell(%(cell_id)s, false);
    \$('#cell_input_%(cell_id)s').val(saved_input);
    send_cell_input(%(cell_id)s);
    cell_input_resize(%(cell_id)s);
""" % locals()

    update_button = r"""<button onclick="
    var f, g, saved_input;
    g = \$('#iframe_graph_editor_%(cell_id)s')[0].contentWindow.update_sage();

    if (g[2] === '') {
        alert('You need to give a Sage variable name to the graph, before saving it.');
        return;
    }
    f = g[2] + ' = Graph(' + g[0] + '); ' + g[2] + '.set_pos(' + g[1] + '); '
    %(eval_strategy)s
">Save</button>""" % locals()

    graph_js = graph_to_js(graph)
    data_fields = """<input type="hidden" id="graph_data_%(cell_id)s" value="%(graph_js)s"><input type="hidden" id="graph_name_%(cell_id)s" value="%(graph_name)s">""" % locals()

    return html(r"""<div id="graph_editor_%(cell_id)s"><table><tbody>
      <tr><td><iframe style="width: 800px; height: 400px; border: 0;" id="iframe_graph_editor_%(cell_id)s" src="/javascript/graph_editor/graph_editor.html?cell_id=%(cell_id)s"></iframe>%(data_fields)s</td></tr>
      <tr><td>%(update_button)s%(close_button)s</td></tr>
</tbody></table></div>""" % locals())

# This is commented out because the mouse_out call raises an error in
# Firebug's console when the event fires but the function itself has
# not yet been loaded.

#      <tr><td><iframe style="width: 800px; height: 400px; border: 0;" id="iframe_graph_editor_%(cell_id)s" src="/javascript/graph_editor/graph_editor.html?cell_id=%(cell_id)s" onmouseout="\$('#iframe_graph_editor_%(cell_id)s')[0].contentWindow.mouse_out();"></iframe>%(data_fields)s</td></tr>
