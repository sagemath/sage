"""
Miscellaneous Functions
"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.all import PermutationGroup, PermutationGroup_generic, PermutationGroupElement, SymmetricGroup
from sage.misc.misc_c import prod
from functools import wraps

def change_support(perm, support, change_perm=None):
    """
    Changes the support of a permutation defined on [1, ..., n] to
    support.

    EXAMPLES::

        sage: from sage.combinat.species.misc import change_support
        sage: p = PermutationGroupElement((1,2,3)); p
        (1,2,3)
        sage: change_support(p, [3,4,5])
        (3,4,5)
    """
    if change_perm is None:
        change_perm = prod([PermutationGroupElement((i+1,support[i])) for i in range(len(support)) if i+1 != support[i]],  PermutationGroupElement([], SymmetricGroup(support)))

    if isinstance(perm, PermutationGroup_generic):
        return PermutationGroup([change_support(g, support, change_perm) for g in perm.gens()])

    return change_perm*perm*~change_perm


def accept_size(f):
    """
    The purpose of this decorator is to change calls like
    species.SetSpecies(size=1) to species.SetSpecies(min=1, max=2).
    This is to make caching species easier and to restrict the number
    of parameters that the lower level code needs to know about.

    EXAMPLES::

        sage: from sage.combinat.species.misc import accept_size
        sage: def f(*args, **kwds):
        ....:       print("{} {}".format(args, sorted(kwds.items())))
        sage: f = accept_size(f)
        sage: f(min=1)
        () [('min', 1)]
        sage: f(size=2)
        () [('max', 3), ('min', 2)]
    """
    @wraps(f)
    def new_func(*args, **kwds):
        if 'size' in kwds:
            if 'min' in kwds or 'max' in kwds:
                raise ValueError("cannot specify both size and (min or max)")
            kwds['min'] = kwds['size']
            kwds['max'] = kwds['size'] + 1
            del kwds['size']
        return f(*args, **kwds)
    return new_func
