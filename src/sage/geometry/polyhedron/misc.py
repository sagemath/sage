# sage_setup: distribution = sagemath-polyhedra
r"""
Miscellaneous helper functions
"""
# **********************************************************************
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# **********************************************************************


def _to_space_separated_string(l, base_ring=None):
    """
    Convert a container to a space-separated string.

    INPUT:

    - ``l`` -- anything iterable

    - ``base_ring`` -- ring (default: ``None``); convert this ring, if given

    OUTPUT: string

    EXAMPLES::

        sage: import sage.geometry.polyhedron.misc as P
        sage: P._to_space_separated_string([2,3])
        '2 3'
        sage: P._to_space_separated_string([2, 1/5], RDF)                               # needs sage.rings.real_double
        '2.0 0.2'
    """
    if base_ring:
        return ' '.join(repr(base_ring(x)) for x in l)
    return ' '.join(repr(x) for x in l)


def _set_to_None_if_empty(x):
    """
    Helper function to clean up arguments.

    This returns None if x is None or x is an empty container.

    EXAMPLES::

        sage: import sage.geometry.polyhedron.misc as P
        sage: None == P._set_to_None_if_empty([])
        True
        sage: P._set_to_None_if_empty([1])
        [1]
    """
    if x is None:
        return x
    x = list(x)
    if not x:
        return None
    return x


def _make_listlist(x):
    """
    Helper function to clean up arguments.

    INPUT:

    - ``x`` -- ``None`` or an iterable of iterables

    OUTPUT: list of lists

    EXAMPLES::

        sage: import sage.geometry.polyhedron.misc as P
        sage: [] == P._make_listlist(tuple())
        True
        sage: [] == P._make_listlist(None)
        True
        sage: P._make_listlist([(1,2),[3,4]])
        [[1, 2], [3, 4]]
    """
    if x is None:
        return []
    return [list(y) for y in x]


def _common_length_of(l1, l2=None, l3=None):
    """
    The arguments are containers or ``None``. The function applies
    ``len()`` to each element, and returns the common length. If the
    length differs, :exc:`ValueError` is raised. Used to check arguments.

    OUTPUT:

    A tuple (number of entries, common length of the entries)

    EXAMPLES::

        sage: import sage.geometry.polyhedron.misc as P
        sage: P._common_length_of([[1,2,3],[1,3,34]])
        (2, 3)
    """
    args = []
    if l1 is not None:
        args.append(l1)
    if l2 is not None:
        args.append(l2)
    if l3 is not None:
        args.append(l3)

    length = None
    num = 0
    for l in args:
        for i in l:
            num += 1
            length_i = len(i)
            if length is not None and length_i != length:
                raise ValueError("Argument lengths differ!")
            length = length_i

    return num, length
