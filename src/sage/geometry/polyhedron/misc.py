r"""
Miscellaneous helper functions
"""


########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################






#########################################################################
def _to_space_separated_string(l):
    """
    Converts a container to a space-separated string.

    INPUT:

    - ``l`` -- anything iterable.

    OUTPUT:

    String.

    EXAMPLES::

        sage: import sage.geometry.polyhedron.misc as P
        sage: P._to_space_separated_string([2,3])
        '2 3'
    """
    s = '';
    for x in l:
        if len(s)>0: s += ' '
        s += repr(x)
    return s


#########################################################################
def _set_to_None_if_empty(x):
    """
    Helper function to clean up arguments: Returns None if x==None or
    x is an empty container.

    EXAMPLES::

        sage: import sage.geometry.polyhedron.misc as P
        sage: None == P._set_to_None_if_empty([])
        True
        sage: P._set_to_None_if_empty([1])
        [1]
    """
    if x is None: return x
    x = list(x)
    if len(x)==0: return None
    return x


#########################################################################
def _make_listlist(x):
    """
    Helper function to clean up arguments.

    INPUT:

    - ``x`` -- ``None`` or an iterable of iterables.

    OUTPUT:

    A list of lists.

    EXAMPLES::

        sage: import sage.geometry.polyhedron.misc as P
        sage: [] == P._make_listlist(tuple())
        True
        sage: [] == P._make_listlist(None)
        True
        sage: P._make_listlist([(1,2),[3,4]])
        [[1, 2], [3, 4]]
    """
    if x is None: return []
    return [list(y) for y in x]


#########################################################################
def _common_length_of(l1, l2=None, l3=None):
    """
    The arguments are containers or ``None``. The function applies
    ``len()`` to each element, and returns the common length. If the
    length differs, ``ValueError`` is raised. Used to check arguments.

    OUTPUT:

    A tuple (number of entries, common length of the entries)

    EXAMPLES::

        sage: import sage.geometry.polyhedron.misc as P
        sage: P._common_length_of([[1,2,3],[1,3,34]])
        (2, 3)
    """
    args = [];
    if l1 is not None: args.append(l1)
    if l2 is not None: args.append(l2)
    if l3 is not None: args.append(l3)

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
