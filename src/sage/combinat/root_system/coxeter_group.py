"""
Coxeter Groups
"""
# ***************************************************************************
#       Copyright (C) 2010 Nicolas Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.combinat.root_system.cartan_type import CartanType
from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.root_system.reflection_group_real', 'ReflectionGroup')
lazy_import('sage.combinat.root_system.weyl_group', 'WeylGroup')


def CoxeterGroup(data, implementation='reflection', base_ring=None, index_set=None):
    """
    Return an implementation of the Coxeter group given by ``data``.

    INPUT:

    - ``data`` -- a Cartan type (or coercible into; see :class:`CartanType`)
      or a Coxeter matrix or graph

    - ``implementation`` -- (default: ``'reflection'``) can be one of
      the following:

      * ``'permutation'`` -- as a permutation representation
      * ``'matrix'`` -- as a Weyl group (as a matrix group acting on the
        root space); if this is not implemented, this uses the "reflection"
        implementation
      * ``'coxeter3'`` -- using the coxeter3 package
      * ``'reflection'`` -- as elements in the reflection representation; see
        :class:`~sage.groups.matrix_gps.coxeter_groups.CoxeterMatrixGroup`

    - ``base_ring`` -- (optional) the base ring for the ``'reflection'``
      implementation

    - ``index_set`` -- (optional) the index set for the ``'reflection'``
      implementation

    EXAMPLES:

    Now assume that ``data`` represents a Cartan type. If
    ``implementation`` is not specified, the reflection representation
    is returned::

        sage: W = CoxeterGroup(["A",2]); W                                              # needs sage.libs.gap
        Finite Coxeter group over Integer Ring with Coxeter matrix:
        [1 3]
        [3 1]

        sage: W = CoxeterGroup(["A",3,1]); W                                            # needs sage.libs.gap
        Coxeter group over Integer Ring with Coxeter matrix:
        [1 3 2 3]
        [3 1 3 2]
        [2 3 1 3]
        [3 2 3 1]

        sage: W = CoxeterGroup(['H',3]); W                                              # needs sage.libs.gap
        Finite Coxeter group over Number Field in a with defining polynomial x^2 - 5
         with a = 2.236067977499790? with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]

    We now use the ``implementation`` option::

        sage: W = CoxeterGroup(["A",2], implementation='permutation'); W    # optional - gap3
        Permutation Group with generators [(1,4)(2,3)(5,6), (1,3)(2,5)(4,6)]
        sage: W.category()                                                  # optional - gap3
        Join of Category of finite enumerated permutation groups
            and Category of finite Weyl groups
            and Category of well generated finite irreducible complex reflection groups

        sage: W = CoxeterGroup(["A",2], implementation='matrix'); W                     # needs sage.libs.gap
        Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)

        sage: W = CoxeterGroup(["H",3], implementation='matrix'); W                     # needs sage.libs.gap sage.rings.number_field
        Finite Coxeter group over Number Field in a with defining polynomial x^2 - 5
         with a = 2.236067977499790? with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]

        sage: W = CoxeterGroup(["H",3], implementation='reflection'); W                 # needs sage.libs.gap sage.rings.number_field
        Finite Coxeter group over Number Field in a with defining polynomial x^2 - 5
         with a = 2.236067977499790? with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]

        sage: W = CoxeterGroup(["A",4,1], implementation='permutation')                 # needs sage.libs.gap
        Traceback (most recent call last):
        ...
        ValueError: the type must be finite

        sage: W = CoxeterGroup(["A",4], implementation='chevie'); W         # optional - gap3
        Irreducible real reflection group of rank 4 and type A4

    We use the different options for the "reflection" implementation::

        sage: W = CoxeterGroup(["H",3], implementation='reflection', base_ring=RR); W   # needs sage.libs.gap
        Finite Coxeter group over Real Field with 53 bits of precision with Coxeter matrix:
        [1 3 2]
        [3 1 5]
        [2 5 1]
        sage: W = CoxeterGroup([[1,10],[10,1]], implementation='reflection',            # needs sage.symbolics
        ....:                  index_set=['a','b'], base_ring=SR); W
        Finite Coxeter group over Symbolic Ring with Coxeter matrix:
        [ 1 10]
        [10  1]

    TESTS::

        sage: W = groups.misc.CoxeterGroup(["H",3])                                     # needs sage.graphs sage.groups
    """
    if implementation not in ["permutation", "matrix", "coxeter3", "reflection", "chevie", None]:
        raise ValueError("invalid type implementation")

    from sage.groups.matrix_gps.coxeter_group import CoxeterMatrixGroup

    try:
        cartan_type = CartanType(data)
    except (TypeError, ValueError): # If it is not a Cartan type, try to see if we can represent it as a matrix group
        return CoxeterMatrixGroup(data, base_ring, index_set)

    if implementation is None:
        implementation = "matrix"

    if implementation == "reflection":
        return CoxeterMatrixGroup(cartan_type, base_ring, index_set)
    if implementation == "coxeter3":
        try:
            from sage.libs.coxeter3.coxeter_group import CoxeterGroup
        except ImportError:
            raise RuntimeError("coxeter3 must be installed")
        else:
            return CoxeterGroup(cartan_type)
    if implementation == "permutation":
        if not cartan_type.is_finite():
            raise ValueError("the type must be finite")
        if cartan_type.is_crystallographic():
            return WeylGroup(cartan_type, implementation='permutation')
        return ReflectionGroup(cartan_type, index_set=index_set)
    elif implementation == "matrix":
        if cartan_type.is_crystallographic():
            return WeylGroup(cartan_type)
        return CoxeterMatrixGroup(cartan_type, base_ring, index_set)
    elif implementation == "chevie":
        return ReflectionGroup(cartan_type, index_set=index_set)

    raise NotImplementedError("Coxeter group of type {} as {} group not implemented".format(cartan_type, implementation))


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.root_system.coxeter_group', 'CoxeterGroupAsPermutationGroup',  ReflectionGroup)
