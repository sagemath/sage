r"""
ReST index of functions

This module contains a function that generates a ReST index table of functions
for use in doc-strings.

{INDEX_OF_FUNCTIONS}

"""

import inspect

from sage.misc.sageinspect import _extract_embedded_position
from sage.misc.sageinspect import is_function_or_cython_function as _isfunction


def gen_rest_table_index(obj, names=None, sort=True, only_local_functions=True):
    r"""
    Return a ReST table describing a list of functions.

    The list of functions can either be given explicitly, or implicitly as the
    functions/methods of a module or class.

    In the case of a class, only non-inherited methods are listed.

    INPUT:

    - ``obj`` -- a list of functions, a module or a class. If given a list of
      functions, the generated table will consist of these. If given a module
      or a class, all functions/methods it defines will be listed, except
      deprecated or those starting with an underscore. In the case of a class,
      note that inherited methods are not displayed.

    - ``names`` -- a dictionary associating a name to a function. Takes
      precedence over the automatically computed name for the functions. Only
      used when ``list_of_entries`` is a list.

    - ``sort`` (boolean; ``True``) -- whether to sort the list of methods
      lexicographically.

    - ``only_local_functions`` (boolean; ``True``) -- if ``list_of_entries`` is
      a module, ``only_local_functions = True`` means that imported functions
      will be filtered out. This can be useful to disable for making indexes of
      e.g. catalog modules such as :mod:`sage.coding.codes_catalog`.

    .. WARNING::

        The ReST tables returned by this function use '@' as a delimiter for
        cells. This can cause trouble if the first sentence in the documentation
        of a function contains the '@' character.

    EXAMPLES::

        sage: from sage.misc.rest_index_of_methods import gen_rest_table_index
        sage: print(gen_rest_table_index([graphs.PetersenGraph]))
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: @
        <BLANKLINE>
           :func:`~sage.graphs.generators.smallgraphs.PetersenGraph` @ Return the Petersen Graph.

    The table of a module::

        sage: print(gen_rest_table_index(sage.misc.rest_index_of_methods))
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: @
        <BLANKLINE>
           :func:`~sage.misc.rest_index_of_methods.doc_index` @ Attribute an index name to a function.
           :func:`~sage.misc.rest_index_of_methods.gen_rest_table_index` @ Return a ReST table describing a list of functions.
           :func:`~sage.misc.rest_index_of_methods.gen_thematic_rest_table_index` @ Return a ReST string of thematically sorted function (or methods) of a module (or class).
           :func:`~sage.misc.rest_index_of_methods.list_of_subfunctions` @ Returns the functions (resp. methods) of a given module (resp. class) with their names.
        <BLANKLINE>
        <BLANKLINE>

    The table of a class::

        sage: print(gen_rest_table_index(Graph))
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: @
        ...
           :meth:`~sage.graphs.graph.Graph.sparse6_string` @ Return the sparse6 representation of the graph as an ASCII string.
        ...

    TESTS:

    When the first sentence of the docstring spans over several lines::

        sage: def a():
        ....:     r'''
        ....:     Here is a very very very long sentence
        ....:     that spans on several lines.
        ....:
        ....:     EXAMP...
        ....:     '''
        ....:     print("hey")
        sage: 'Here is a very very very long sentence that spans on several lines' in gen_rest_table_index([a])
        True

    The inherited methods do not show up::

        sage: gen_rest_table_index(sage.combinat.posets.lattices.FiniteLatticePoset).count('\n') < 75
        True
        sage: from sage.graphs.generic_graph import GenericGraph
        sage: A = gen_rest_table_index(Graph).count('\n')
        sage: B = gen_rest_table_index(GenericGraph).count('\n')
        sage: A < B
        True

    When ``only_local_functions`` is ``False``, we do not include
    ``gen_rest_table_index`` itself::

        sage: print(gen_rest_table_index(sage.misc.rest_index_of_methods, only_local_functions=True))
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: @
        <BLANKLINE>
           :func:`~sage.misc.rest_index_of_methods.doc_index` @ Attribute an index name to a function.
           :func:`~sage.misc.rest_index_of_methods.gen_rest_table_index` @ Return a ReST table describing a list of functions.
           :func:`~sage.misc.rest_index_of_methods.gen_thematic_rest_table_index` @ Return a ReST string of thematically sorted function (or methods) of a module (or class).
           :func:`~sage.misc.rest_index_of_methods.list_of_subfunctions` @ Returns the functions (resp. methods) of a given module (resp. class) with their names.
        <BLANKLINE>
        <BLANKLINE>
        sage: print(gen_rest_table_index(sage.misc.rest_index_of_methods, only_local_functions=False))
        .. csv-table::
           :class: contentstable
           :widths: 30, 70
           :delim: @
        <BLANKLINE>
           :func:`~sage.misc.rest_index_of_methods.doc_index` @ Attribute an index name to a function.
           :func:`~sage.misc.rest_index_of_methods.gen_thematic_rest_table_index` @ Return a ReST string of thematically sorted function (or methods) of a module (or class).
           :func:`~sage.misc.rest_index_of_methods.list_of_subfunctions` @ Returns the functions (resp. methods) of a given module (resp. class) with their names.
        <BLANKLINE>
        <BLANKLINE>

    A function that is imported into a class under a different name is listed
    under its 'new' name::

        sage: 'cliques_maximum' in gen_rest_table_index(Graph)
        True
        sage: 'all_max_cliques`' in gen_rest_table_index(Graph)
        False
    """
    if names is None:
        names = {}

    # If input is a class/module, we list all its non-private and methods/functions
    if inspect.isclass(obj) or inspect.ismodule(obj):
        list_of_entries, names = list_of_subfunctions(
                obj, only_local_functions=only_local_functions)
    else:
        list_of_entries = obj

    fname = lambda x: names.get(x, getattr(x, "__name__", ""))

    assert isinstance(list_of_entries, list)

    s = [".. csv-table::",
         "   :class: contentstable",
         "   :widths: 30, 70",
         "   :delim: @\n"]

    if sort:
        list_of_entries.sort(key=fname)

    for e in list_of_entries:
        if inspect.ismethod(e):
            link = ":meth:`~{module}.{cls}.{func}`".format(
                module=e.im_class.__module__, cls=e.im_class.__name__,
                func=fname(e))
        elif _isfunction(e) and inspect.isclass(obj):
            link = ":meth:`~{module}.{cls}.{func}`".format(
                module=obj.__module__, cls=obj.__name__, func=fname(e))
        elif _isfunction(e):
            link = ":func:`~{module}.{func}`".format(
                module=e.__module__, func=fname(e))
        else:
            continue

        # Extract lines injected by cython
        doc = e.__doc__
        doc_tmp = _extract_embedded_position(doc)
        if doc_tmp:
            doc = doc_tmp[0]

        # Descriptions of the method/function
        if doc:
            desc = doc.split('\n\n')[0]                             # first paragraph
            desc = " ".join(x.strip() for x in desc.splitlines())   # concatenate lines
            desc = desc.strip()                                     # remove leading spaces
        else:
            desc = "NO DOCSTRING"

        s.append("   {} @ {}".format(link, desc.lstrip()))

    return '\n'.join(s) + '\n'


def list_of_subfunctions(root, only_local_functions=True):
    r"""
    Returns the functions (resp. methods) of a given module (resp. class) with their names.

    INPUT:

    - ``root`` -- the module, or class, whose elements are to be listed.

    - ``only_local_functions`` (boolean; ``True``) -- if ``root`` is a module,
      ``only_local_functions = True`` means that imported functions will be
      filtered out. This can be useful to disable for making indexes of
      e.g. catalog modules such as :mod:`sage.coding.codes_catalog`.

    OUTPUT:

    A pair ``(list,dict)`` where ``list`` is a list of function/methods and
    ``dict`` associates to every function/method the name under which it appears
    in ``root``.

    EXAMPLES::

        sage: from sage.misc.rest_index_of_methods import list_of_subfunctions
        sage: l = list_of_subfunctions(Graph)[0]
        sage: Graph.bipartite_color in l
        True

    TESTS:

    A ``staticmethod`` is not callable. We must handle them correctly, however::

        sage: class A:
        ....:     x = staticmethod(Graph.order)
        sage: list_of_subfunctions(A)
        ([<function GenericGraph.order at 0x...>],
         {<function GenericGraph.order at 0x...>: 'x'})

    """
    if inspect.ismodule(root):
        ismodule = True
    elif inspect.isclass(root):
        ismodule = False
        superclasses = inspect.getmro(root)[1:]
    else:
        raise ValueError("'root' must be a module or a class.")

    def local_filter(f,name):
        if only_local_functions:
            if ismodule:
                return inspect.getmodule(root) == inspect.getmodule(f)
            else:
                return not any(hasattr(s,name) for s in superclasses)
        else:
            return inspect.isclass(root) or not (f is gen_rest_table_index)

    functions =  {getattr(root,name):name for name,f in root.__dict__.items() if
                  (not name.startswith('_')          and # private functions
                   not hasattr(f,'trac_number')      and # deprecated functions
                   not inspect.isclass(f)            and # classes
                   callable(getattr(f,'__func__',f)) and # e.g. GenericGraph.graphics_array_defaults
                   local_filter(f,name))                 # possibly filter imported functions
                  }

    return list(functions.keys()), functions


def gen_thematic_rest_table_index(root,additional_categories=None,only_local_functions=True):
    r"""
    Return a ReST string of thematically sorted function (or methods) of a module (or class).

    INPUT:

    - ``root`` -- the module, or class, whose elements are to be listed.

    - ``additional_categories`` -- a dictionary associating a category (given as
      a string) to a function's name. Can be used when the decorator
      :func:`doc_index` does not work on a function.

    - ``only_local_functions`` (boolean; ``True``) -- if ``root`` is a module,
      ``only_local_functions = True`` means that imported functions will be
      filtered out. This can be useful to disable for making indexes of
      e.g. catalog modules such as :mod:`sage.coding.codes_catalog`.

    EXAMPLES::

        sage: from sage.misc.rest_index_of_methods import gen_thematic_rest_table_index, list_of_subfunctions
        sage: l = list_of_subfunctions(Graph)[0]
        sage: Graph.bipartite_color in l
        True
    """
    from collections import defaultdict
    if additional_categories is None:
        additional_categories = {}

    functions, names = list_of_subfunctions(root,
                                            only_local_functions=only_local_functions)
    theme_to_function = defaultdict(list)
    for f in functions:
        if hasattr(f, 'doc_index'):
            doc_ind = f.doc_index
        else:
            try:
                doc_ind = additional_categories.get(f.__name__,
                                                    "Unsorted")
            except AttributeError:
                doc_ind = "Unsorted"
        theme_to_function[doc_ind].append(f)
    s = ["**"+theme+"**\n\n"+gen_rest_table_index(list_of_functions,names=names)
         for theme, list_of_functions in sorted(theme_to_function.items())]
    return "\n\n".join(s)


def doc_index(name):
    r"""
    Attribute an index name to a function.

    This decorator can be applied to a function/method in order to specify in
    which index it must appear, in the index generated by
    :func:`gen_thematic_rest_table_index`.

    INPUT:

    - ``name`` -- a string, which will become the title of the index in which
      this function/method will appear.

    EXAMPLES::

        sage: from sage.misc.rest_index_of_methods import doc_index
        sage: @doc_index("Wouhouuuuu")
        ....: def a():
        ....:     print("Hey")
        sage: a.doc_index
        'Wouhouuuuu'
    """
    def hey(f):
        setattr(f,"doc_index",name)
        return f
    return hey

__doc__ = __doc__.format(INDEX_OF_FUNCTIONS=gen_rest_table_index([gen_rest_table_index]))
