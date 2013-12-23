r"""
Class inheritance graphs
"""
#*****************************************************************************
#  Copyright (C) 2007 William Stein <wstein@math.ucsd.edu>
#                2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

import inspect

def class_graph(top, depth=5, name_filter=None, classes=None, as_graph = True):
    """
    Returns the class inheritance graph of a module, class, or object

    INPUT:

     - ``top`` -- the module, class, or object to start with (e.g. ``sage``, ``Integer``, ``3``)
     - ``depth`` -- maximal recursion depth within submodules (default: 5)
     - ``name_filter`` -- e.g. 'sage.rings' to only consider classes in :mod:`sage.rings`
     - ``classes`` -- optional dictionary to be filled in (it is also returned)
     - ``as_graph`` -- a boolean (default: True)

    OUTPUT:

     - An oriented graph, with class names as vertices, and an edge
       from each class to each of its bases.

    EXAMPLES:

    We construct the inheritance graph of the classes within a given module::

        sage: from sage.rings.polynomial.padics import polynomial_padic_capped_relative_dense, polynomial_padic_flat
        sage: G = class_graph(sage.rings.polynomial.padics); G
        Digraph on 6 vertices
        sage: G.vertices()
        ['Polynomial', 'Polynomial_generic_dense', 'Polynomial_generic_domain', 'Polynomial_padic', 'Polynomial_padic_capped_relative_dense', 'Polynomial_padic_flat']
        sage: G.edges(labels=False)
        [('Polynomial_padic', 'Polynomial'), ('Polynomial_padic_capped_relative_dense', 'Polynomial_generic_domain'), ('Polynomial_padic_capped_relative_dense', 'Polynomial_padic'), ('Polynomial_padic_flat', 'Polynomial_generic_dense'), ('Polynomial_padic_flat', 'Polynomial_padic')]

    We construct the inheritance graph of a given class::

        sage: class_graph(Parent).edges(labels=False)
        [('CategoryObject', 'SageObject'), ('Parent', 'CategoryObject'), ('SageObject', 'object')]

    We construct the inheritance graph of the class of an object::

        sage: class_graph([1,2,3]).edges(labels=False)
        [('list', 'object')]

    .. warning:: the output of ``class_graph`` used to be a dictionary
       mapping each class name to the list of names of its bases. This
       can be emulated by setting the option ``as_graph`` to ``False``::

        sage: class_graph(sage.rings.polynomial.padics, depth=2, as_graph=False)
        {'Polynomial_padic': ['Polynomial'],
        'Polynomial_padic_capped_relative_dense': ['Polynomial_generic_domain', 'Polynomial_padic'],
        'Polynomial_padic_flat': ['Polynomial_generic_dense', 'Polynomial_padic']}

    .. note:: the ``classes`` and ``as_graph`` options are mostly
       intended for internal recursive use.

    .. note:: ``class_graph`` does not yet handle nested classes

    TESTS::

        sage: G = class_graph(sage.rings.polynomial.padics, depth=2); G
        Digraph on 6 vertices
    """
    # This function descends recursively down the submodules of the
    # top module (if ``top`` is a module) and then down the hierarchy
    # of classes. Along the way, the result is stored in the "global"
    # dictionary ``classes`` which associates to each class the list
    # of its bases.

    # Termination
    if depth < 0:
        return classes

    # (first recursive call)
    if classes is None:
        classes = dict()

    # Build the list ``children`` of submodules (resp. base classes)
    # of ``top`` the function will recurse through
    if inspect.ismodule(top):
        if top.__name__.endswith('.all'): # Ignore sage.rings.all and friends
            return classes
        if name_filter is None:
            name_filter = top.__name__
        children = [item for item in top.__dict__.values()
                       if inspect.ismodule(item) or inspect.isclass(item)]
        depth = depth -1
    elif inspect.isclass(top):
        if name_filter is None:
            name_filter = ""
        if not top.__module__.startswith(name_filter):
            return classes
        children = top.__bases__
        classes[top.__name__] = [e.__name__ for e in children]
    else: # top is a plain Python object; inspect its class
        children = [top.__class__]

    # Recurse
    for child in children:
        class_graph(child, depth = depth, name_filter=name_filter, classes=classes, as_graph = False)

    # (first recursive call): construct the graph
    if as_graph:
        from sage.graphs.digraph import DiGraph
        return DiGraph(classes)
    else:
        return classes

