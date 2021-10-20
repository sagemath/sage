r"""
`\ZZ`-Filtered Vector Spaces

This module implements filtered vector spaces, that is, a descending
sequence of vector spaces

.. MATH::

    \cdots \supset F_d \supset F_{d+1} \supset F_{d+2} \supset \cdots

with degrees `d\in \ZZ`. It is not required that `F_d` is the entire
ambient space for `d\ll 0` (see
:meth:`~FilteredVectorSpace_class.is_exhaustive`) nor that `F_d=0` for
`d\gg 0` (see :meth:`~FilteredVectorSpace_class.is_separating`). To
construct a filtered vector space, use the :func:`FilteredVectorSpace`
command. It supports easy creation of simple filtrations, for example
the trivial one::

    sage: FilteredVectorSpace(2, base_ring=RDF)
    RDF^2

The next-simplest filtration has a single non-trivial inclusion
between `V_d` and `V_{d+1}`::

    sage: d = 1
    sage: V = FilteredVectorSpace(2, d);  V
    QQ^2 >= 0
    sage: [V.get_degree(i).dimension() for i in range(0,4)]
    [2, 2, 0, 0]

To construct general filtrations, you need to tell Sage about generating
vectors for the nested subspaces. For example, a dictionary whose keys
are the degrees and values are a list of generators::

    sage: r1 = (1, 0, 5)
    sage: r2 = (0, 1, 2)
    sage: r3 = (1, 2, 1)
    sage: V = FilteredVectorSpace({0:[r1, r2, r3], 1:[r1, r2], 3:[r1]});  V
    QQ^3 >= QQ^2 >= QQ^1 >= QQ^1 >= 0

For degrees `d` that are not specified, the associated vector subspace
is the same as the next-lower degree, that is, `V_d \simeq
V_{d-1}`. In the above example, this means that

* `V_d \simeq \QQ^3` for `d<0`
* `V_0 = \mathop{span}(r_1, r_2) \simeq \QQ^2`
* `V_1 = V_2 = \mathop{span}(r_3) \simeq \QQ`
* `V_d = 0` for `d \geq 3`

That is::

    sage: V.get_degree(0) == V
    True
    sage: V.get_degree(1) == V.span([r1, r2])
    True
    sage: V.get_degree(2) == V.get_degree(3) == V.span([r1])
    True
    sage: V.get_degree(4) == V.get_degree(5) == V.span([])
    True

If you have many generators you can just pass the generators once and
then refer to them by index::

    sage: FilteredVectorSpace([r1, r2, r3], {0:[0,1,2], 1:[1,2], 3:[1]})
    QQ^3 >= QQ^2 >= QQ^1 >= QQ^1 >= 0

Note that generators for the degree-`d` subspace of the filtration are
automatically generators for all lower degrees. For example, here we
do not have to specify the ray `r_2` separately in degree 1::

    sage: FilteredVectorSpace([r1, r2, r3], {0:[0   ], 1:[1]})
    QQ^2 >= QQ^1 >= 0 in QQ^3
    sage: FilteredVectorSpace([r1, r2, r3], {0:[0, 1], 1:[1]})
    QQ^2 >= QQ^1 >= 0 in QQ^3

The degree can be infinite (plus infinity), this allows construction
of filtered vector spaces that are not eventually zero in high
degree::

    sage: FilteredVectorSpace([r1, r2, r3], {0:[0,1], oo:[1]})
    QQ^2 >= QQ^1 in QQ^3

Any field can be used as the vector space base. For example a finite
field::

    sage: F.<a> = GF(5^3)
    sage: r1 = (a, 0, F(5));  r1
    (a, 0, 0)
    sage: FilteredVectorSpace([r1, r2, r3], {0:[0,1], oo:[1]}, base_ring=F)
    GF(125)^2 >= GF(125)^1 in GF(125)^3

Or the algebraic field::

    sage: r1 = (1, 0, 1+QQbar(I));  r1
    (1, 0, I + 1)
    sage: FilteredVectorSpace([r1, r2, r3], {0:[0,1], oo:[1]}, base_ring=QQbar)
    Vector space of dimension 2 over Algebraic Field
    >= Vector space of dimension 1 over Algebraic Field
    in Vector space of dimension 3 over Algebraic Field
"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import QQ, ZZ, RDF, RR, Integer
from sage.rings.infinity import InfinityRing, infinity, minus_infinity
from sage.categories.fields import Fields
from sage.modules.free_module import FreeModule_ambient_field, VectorSpace
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method


def is_FilteredVectorSpace(X):
    """
    Test whether ``X`` is a filtered vector space.

    This function is for library use only.

    INPUT:

    - ``X`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: from sage.modules.filtered_vector_space import is_FilteredVectorSpace
        sage: V = FilteredVectorSpace(2, 1)
        sage: is_FilteredVectorSpace(V)
        True
        sage: is_FilteredVectorSpace('ceci n\'est pas une pipe')
        False
    """
    return isinstance(X, FilteredVectorSpace_class)


def FilteredVectorSpace(arg1, arg2=None, base_ring=QQ, check=True):
    r"""
    Construct a filtered vector space.

    INPUT:

    This function accepts various input that determines the vector space and filtration.

    - Just the dimensionFilteredVectorSpace(dimension): Return the trivial filtration
      (where all vector spaces are isomorphic).

    - Dimension and maximal degree, see
      :func:`constructor_from_dim_degree` for arguments. Construct a
      filtration with only one non-trivial step `V\supset 0` at the
      given cutoff degree.

    - A dictionary containing the degrees as keys and a list of vector
      space generators as values, see
      :func:`FilteredVectorSpace_from_generators`

    - Generators and a dictionary containing the degrees as keys and
      the indices of vector space generators as values, see
      :func:`FilteredVectorSpace_from_generators_indices`

    In addition, the following keyword arguments are supported:

    - ``base_ring`` -- a field (optional, default `\QQ`). The base
      field of the vector space. Must be a field.

    EXAMPLES:

    Just the dimension for the trivial filtration::

        sage: FilteredVectorSpace(2)
        QQ^2

    Dimension and degree::

        sage: FilteredVectorSpace(2, 1)
        QQ^2 >= 0

    Dictionary of generators::

        sage: FilteredVectorSpace({1:[(1,0), (0,1)], 3:[(1,0)]})
        QQ^2 >= QQ^1 >= QQ^1 >= 0

    Generators and a dictionary referring to them by index::

        sage: FilteredVectorSpace([(1,0), (0,1)], {1:[0,1], 3:[0]})
        QQ^2 >= QQ^1 >= QQ^1 >= 0
    """
    if base_ring not in Fields():
        raise ValueError('the base_ring argument must be a field')
    if arg1 in ZZ:
        return construct_from_dim_degree(arg1, arg2, base_ring, check)
    elif arg2 is None:
        return construct_from_generators(arg1, base_ring, check)
    else:
        return construct_from_generators_indices(arg1, arg2, base_ring, check)


def normalize_degree(deg):
    """
    Normalized the degree

    - ``deg`` -- something that defines the degree (either integer or
      infinity).

    OUTPUT:

    Plus/minus infinity or a Sage integer.

    EXAMPLES::

        sage: from sage.modules.filtered_vector_space import normalize_degree
        sage: type(normalize_degree(int(1)))
        <class 'sage.rings.integer.Integer'>
        sage: normalize_degree(oo)
        +Infinity
    """
    try:
        return ZZ(deg)
    except TypeError:
        pass
    deg = InfinityRing(deg)
    if deg == infinity:
        return infinity
    if deg == minus_infinity:
        return minus_infinity
    raise ValueError('not integer or infinity')


def construct_from_dim_degree(dim, max_degree, base_ring, check):
    """
    Construct a filtered vector space.

    INPUT:

    - ``dim`` -- integer. The dimension.

    - ``max_degree`` -- integer or infinity. The maximal degree where
      the vector subspace of the filtration is still the entire space.

    EXAMPLES::

        sage: V = FilteredVectorSpace(2, 5);  V
        QQ^2 >= 0
        sage: V.get_degree(5)
        Vector space of degree 2 and dimension 2 over Rational Field
        Basis matrix:
        [1 0]
        [0 1]
        sage: V.get_degree(6)
        Vector space of degree 2 and dimension 0 over Rational Field
        Basis matrix:
        []

        sage: FilteredVectorSpace(2, oo)
        QQ^2
        sage: FilteredVectorSpace(2, -oo)
        0 in QQ^2

    TESTS::

        sage: from sage.modules.filtered_vector_space import construct_from_dim_degree
        sage: V = construct_from_dim_degree(2, 5, QQ, True);  V
        QQ^2 >= 0
    """
    if dim not in ZZ:
        raise ValueError('dimension must be an integer')
    dim = ZZ(dim)
    from sage.matrix.constructor import identity_matrix
    generators = identity_matrix(base_ring, dim).columns()
    filtration = dict()
    if max_degree is None:
        max_degree = infinity
    filtration[normalize_degree(max_degree)] = range(dim)
    return construct_from_generators_indices(generators, filtration, base_ring, check)


def construct_from_generators(filtration, base_ring, check):
    """
    Construct a filtered vector space.

    INPUT:

    - ``filtration`` -- a dictionary of filtration steps. Each
      filtration step is a pair consisting of an integer degree and a
      list/tuple/iterable of vector space generators. The integer
      ``degree`` stipulates that all filtration steps of degree higher
      or equal than ``degree`` (up to the next filtration step) are
      said subspace.

    EXAMPLES::

        sage: from sage.modules.filtered_vector_space import construct_from_generators
        sage: r = [1, 2]
        sage: construct_from_generators({1:[r]}, QQ, True)
        QQ^1 >= 0 in QQ^2
    """
    def normalize_gen(v):
        return tuple(map(base_ring, v))

    # convert generator notation to generator+indices
    if len(filtration) == 0:
        raise ValueError('you need to specify at least one ray to deduce the dimension')
    generators = set()
    for gens in filtration.values():
        generators.update(normalize_gen(g) for g in gens)
    generators = tuple(sorted(generators))

    # normalize filtration data
    normalized = dict()
    for deg, gens_deg in filtration.items():
        indices = [generators.index(normalize_gen(v)) for v in gens_deg]
        normalized[deg] = tuple(indices)
    return construct_from_generators_indices(generators, normalized, base_ring, check)


def construct_from_generators_indices(generators, filtration, base_ring, check):
    """
    Construct a filtered vector space.

    INPUT:

    - ``generators`` -- a list/tuple/iterable of vectors, or something
      convertible to them. The generators spanning various
      subspaces.

    - ``filtration`` -- a list or iterable of filtration steps. Each
      filtration step is a pair ``(degree, ray_indices)``. The
      ``ray_indices`` are a list or iterable of ray indices, which
      span a subspace of the vector space. The integer ``degree``
      stipulates that all filtration steps of degree higher or equal
      than ``degree`` (up to the next filtration step) are said
      subspace.

    EXAMPLES::

        sage: from sage.modules.filtered_vector_space import construct_from_generators_indices
        sage: gens = [(1,0), (0,1), (-1,-1)]
        sage: V = construct_from_generators_indices(gens, {1:[0,1], 3:[1]}, QQ, True);  V
        QQ^2 >= QQ^1 >= QQ^1 >= 0

    TESTS::

        sage: gens = [(int(1),int(0)), (0,1), (-1,-1)]
        sage: construct_from_generators_indices(iter(gens), {int(0):[0, int(1)], 2:[2]}, QQ, True)
        QQ^2 >= QQ^1 >= QQ^1 >= 0
    """
    # normalize generators
    generators = [list(g) for g in generators]

    # deduce dimension
    if len(generators) == 0:
        dim = ZZ(0)
    else:
        dim = ZZ(len(generators[0]))
    ambient = VectorSpace(base_ring, dim)

    # complete generators to a generating set
    if matrix(base_ring, generators).rank() < dim:
        complement = ambient.span(generators).complement()
        generators = generators + list(complement.gens())
    # normalize generators II
    generators = tuple(ambient(v) for v in generators)

    for v in generators:
        v.set_immutable()

    # normalize filtration data
    normalized = dict()
    for deg, gens in filtration.items():
        deg = normalize_degree(deg)
        gens = [ZZ(i) for i in gens]
        if any(i < 0 or i >= len(generators) for i in gens):
            raise ValueError('generator index out of bounds')
        normalized[deg] = tuple(sorted(gens))
    try:
        del normalized[minus_infinity]
    except KeyError:
        pass
    filtration = normalized

    return FilteredVectorSpace_class(base_ring, dim, generators, filtration, check=check)




class FilteredVectorSpace_class(FreeModule_ambient_field):

    def __init__(self, base_ring, dim, generators, filtration, check=True):
        r"""
        A descending filtration of a vector space

        INPUT:

        - ``base_ring`` -- a field. The base field of the ambient vector space.

        - ``dim`` -- integer. The dimension of the ambient vector space.

        - ``generators`` -- tuple of generators for the ambient vector
          space. These will be used to span the subspaces of the
          filtration.

        - ``filtration`` -- a dictionary of filtration steps in ray
          index notation. See
          :func:`construct_from_generators_indices` for details.

        - ``check`` -- boolean (optional; default: ``True``). Whether
          to perform consistency checks.

        TESTS::

            sage: from sage.modules.filtered_vector_space import FilteredVectorSpace_class
            sage: gens = [(1,0,0), (1,1,0), (1,2,0), (-1,-1, 0), (0,0,1)]
            sage: FilteredVectorSpace_class(QQ, 3, gens, {2:(0,1), oo:(4,)})
            QQ^3 >= QQ^1
            sage: FilteredVectorSpace_class(QQ, 3, gens, {2:(0,1), 3:(4,)})
            QQ^3 >= QQ^1 >= 0

        The trivial filtration::

            sage: FilteredVectorSpace_class(QQ, 3, gens, {}, QQ)
            0 in QQ^3

        The empty vector space::

            sage: FilteredVectorSpace_class(QQ, 0, [], {})
            0

        Higher-degree generators are automatically generators in lower degrees::

            sage: FilteredVectorSpace_class(QQ, 3, gens, {2:(4,), 3:(1,)})
            QQ^2 >= QQ^1 >= 0 in QQ^3
        """
        if check:
            assert isinstance(dim, Integer)
            assert base_ring in Fields()
        super(FilteredVectorSpace_class, self).__init__(base_ring, dim)

        if check:
            assert matrix(generators).rank() == self.dimension()
            assert isinstance(filtration, dict)
            for degree, indices in filtration.items():
                assert isinstance(degree, Integer) or degree == infinity
                assert isinstance(indices, tuple)
                assert all(isinstance(r, Integer) for r in indices)

        # Construct subspaces from the generators and store in self._filt
        def make_subspace(indices):
            return self.span([generators[i] for i in indices])

        indices = set(filtration.pop(infinity, []))
        V = make_subspace(indices)
        filtered_subspaces = [(infinity, V)]
        for deg in reversed(sorted(filtration.keys())):
            next_V = V
            indices.update(filtration[deg])
            V = make_subspace(indices)
            if V == next_V:   # skip trivial filtrations
                continue
            filtered_subspaces.append((deg, V))
        filtered_subspaces.append((minus_infinity, V))
        filtered_subspaces.reverse()
        self._filt = tuple(filtered_subspaces)
        assert self._filt[0][0] is minus_infinity

    def change_ring(self, base_ring):
        """
        Return the same filtration over a different base ring.

        INPUT:

        - ``base_ring`` -- a ring. The new base ring.

        OUTPUT:

        This method returns a new filtered vector space whose
        subspaces are defined by the same generators but over a
        different base ring.

        EXAMPLES::

            sage: V = FilteredVectorSpace(1, 0);  V
            QQ^1 >= 0
            sage: V.change_ring(RDF)
            RDF^1 >= 0
        """
        generators, filtration = self.presentation()
        return FilteredVectorSpace(generators, filtration, base_ring=base_ring)

    def ambient_vector_space(self):
        """
        Return the ambient (unfiltered) vector space.

        OUTPUT:

        A vector space.

        EXAMPLES::

            sage: V = FilteredVectorSpace(1, 0)
            sage: V.ambient_vector_space()
            Vector space of dimension 1 over Rational Field
        """
        return VectorSpace(self.base_ring(), self.dimension())

    @cached_method
    def is_constant(self):
        """
        Return whether the filtration is constant.

        OUTPUT:

        Boolean. Whether the filtered vector spaces are identical in
        all degrees.

        EXAMPLES::

            sage: V = FilteredVectorSpace(2); V
            QQ^2
            sage: V.is_constant()
            True

            sage: V = FilteredVectorSpace(1, 0);  V
            QQ^1 >= 0
            sage: V.is_constant()
            False

            sage: V = FilteredVectorSpace({0:[(1,)]});  V
            QQ^1 >= 0
            sage: V.is_constant()
            False
        """
        f = self._filt
        return (len(f) == 1) or (len(f) == 2 and f[1][0] == infinity)

    def is_exhaustive(self):
        r"""
        Return whether the filtration is exhaustive.

        A filtration $\{F_d\}$ in an ambient vector space $V$ is
        exhaustive if $\cup F_d = V$. See also :meth:`is_separating`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: F = FilteredVectorSpace({0:[(1,1)]});  F
            QQ^1 >= 0 in QQ^2
            sage: F.is_exhaustive()
            False
            sage: G = FilteredVectorSpace(2, 0);  G
            QQ^2 >= 0
            sage: G.is_exhaustive()
            True
        """
        return self.get_degree(minus_infinity).dimension() == \
            self.ambient_vector_space().dimension()

    def is_separating(self):
        r"""
        Return whether the filtration is separating.

        A filtration $\{F_d\}$ in an ambient vector space $V$ is
        exhaustive if $\cap F_d = 0$. See also :meth:`is_exhaustive`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: F = FilteredVectorSpace({0:[(1,1)]});  F
            QQ^1 >= 0 in QQ^2
            sage: F.is_separating()
            True
            sage: G = FilteredVectorSpace({0:[(1,1,0)], oo:[(0,0,1)]});  G
            QQ^2 >= QQ^1 in QQ^3
            sage: G.is_separating()
            False
        """
        return self.get_degree(infinity).dimension() == 0

    @cached_method
    def support(self):
        """
        Return the degrees in which there are non-trivial generators.

        OUTPUT:

        A tuple of integers (and plus infinity) in ascending
        order. The last entry is plus infinity if and only if the
        filtration is not separating (see :meth:`is_separating`).

        EXAMPLES::

            sage: G = FilteredVectorSpace({0:[(1,1,0)], 3:[(0,1,0)]});  G
            QQ^2 >= QQ^1 >= QQ^1 >= QQ^1 >= 0 in QQ^3
            sage: G.support()
            (0, 3)

            sage: G = FilteredVectorSpace({0:[(1,1,0)], 3:[(0,1,0)], oo:[(0,0,1)]});  G
            QQ^3 >= QQ^2 >= QQ^2 >= QQ^2 >= QQ^1
            sage: G.support()
            (0, 3, +Infinity)
        """
        if self.is_separating():
            filt = self._filt[1:-1]
        else:
            filt = self._filt[1:]
        return tuple(f[0] for f in filt)

    @cached_method
    def min_degree(self):
        r"""
        Return the lowest degree of the filtration.

        OUTPUT:

        Integer or plus infinity. The largest degree `d` of the
        (descending) filtration such that the filtered vector space
        `F_d` is still equal to `F_{-\infty}`.

        EXAMPLES::

            sage: FilteredVectorSpace(1, 3).min_degree()
            3
            sage: FilteredVectorSpace(2).min_degree()
            +Infinity
        """
        if self.is_constant():
            return infinity
        return self._filt[1][0]

    @cached_method
    def max_degree(self):
        r"""
        Return the highest degree of the filtration.

        OUTPUT:

        Integer or minus infinity. The smallest degree of the
        filtration such that the filtration is constant to the right.

        EXAMPLES::

            sage: FilteredVectorSpace(1, 3).max_degree()
            4
            sage: FilteredVectorSpace({0:[[1]]}).max_degree()
            1
            sage: FilteredVectorSpace(3).max_degree()
            -Infinity
        """
        f = self._filt
        if len(f) == 1:
            return minus_infinity
        d = f[-1][0]
        if d == infinity:
            if len(f) == 1:
                return minus_infinity
            else:
                return f[-2][0] + 1
        else:
            return d + 1

    def get_degree(self, d):
        r"""
        Return the degree-``d`` entry of the filtration.

        INPUT:

        - ``d`` -- Integer. The desired degree of the filtration.

        OUTPUT:

        The degree-``d`` vector space in the filtration as subspace of
        the ambient space.

        EXAMPLES::

            sage: rays = [(1,0), (1,1), (1,2), (-1,-1)]
            sage: F = FilteredVectorSpace(rays, {3:[1], 1:[1,2]})
            sage: F.get_degree(2)
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 1]
            sage: F.get_degree(oo)
            Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: F.get_degree(-oo)
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
        """
        d = normalize_degree(d)
        for deg, Vdeg in self._filt:
            if d <= deg:
                return Vdeg
        assert False  # unreachable

    def graded(self, d):
        r"""
        Return the associated graded vectorspace.

        INPUT:

        - ``d`` -- integer. The degree.

        OUTPUT:

        The quotient `G_d = F_d / F_{d+1}`.

        EXAMPLES::

            sage: rays = [(1,0), (1,1), (1,2)]
            sage: F = FilteredVectorSpace(rays, {3:[1], 1:[1,2]})
            sage: F.graded(1)
            Vector space quotient V/W of dimension 1 over Rational Field where
            V: Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            W: Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 1]
        """
        return self.get_degree(d).quotient(self.get_degree(d+1))

    def presentation(self):
        """
        Return a presentation in term of generators of various degrees.

        OUTPUT:

        A pair consisting of generators and a filtration suitable as
        input to :func:`~construct_from_generators_indices`.

        EXAMPLES::

            sage: rays = [(1,0), (1,1), (1,2), (-1,-1)]
            sage: F = FilteredVectorSpace(rays, {0:[1, 2], 2:[3]});  F
            QQ^2 >= QQ^1 >= QQ^1 >= 0
            sage: F.presentation()
            (((0, 1), (1, 0), (1, 1)), {0: (1, 0), 2: (2,), +Infinity: ()})
        """
        # this could be done more efficiently with (potentially) less generators
        generators = set()
        filt = self._filt[1:]
        for d, V in filt:
            generators.update(V.echelonized_basis())
        generators = tuple(sorted(generators))

        filtration = dict()
        for d, V in filt:
            indices = [ZZ(generators.index(v)) for v in V.echelonized_basis()]
            filtration[d] = tuple(indices)
        return generators, filtration

    def _repr_field_name(self):
        """
        Return an abbreviated field name as string

        RAISES:

        ``NotImplementedError``: The field does not have an
        abbreviated name defined.

        EXAMPLES::

            sage: FilteredVectorSpace(2, base_ring=QQ)._repr_field_name()
            'QQ'

            sage: F.<a> = GF(9)
            sage: FilteredVectorSpace(2, base_ring=F)._repr_field_name()
            'GF(9)'

            sage: FilteredVectorSpace(2, base_ring=AA)._repr_field_name()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if self.base_ring() == QQ:
            return 'QQ'
        elif self.base_ring() == RDF:
            return 'RDF'
        elif self.base_ring() == RR:
            return 'RR'
        from sage.categories.finite_fields import FiniteFields
        if self.base_ring() in FiniteFields():
            return 'GF({0})'.format(len(self.base_ring()))
        else:
            raise NotImplementedError()

    def _repr_vector_space(self, dim):
        """
        Return a string representation of the vector space of given dimension

        INPUT:

        - ``dim`` -- integer.

        OUTPUT:

        String representation of the vector space of dimension ``dim``.

        EXAMPLES::

            sage: F = FilteredVectorSpace(3, base_ring=RDF)
            sage: F._repr_vector_space(1234)
            'RDF^1234'
            sage: F3 = FilteredVectorSpace(3, base_ring=GF(3))
            sage: F3._repr_vector_space(1234)
            'GF(3)^1234'
            sage: F3 = FilteredVectorSpace(3, base_ring=AA)
            sage: F3._repr_vector_space(1234)
            'Vector space of dimension 1234 over Algebraic Real Field'
        """
        if dim == 0:
            return '0'
        try:
            return self._repr_field_name() + '^' + str(dim)
        except NotImplementedError:
            return repr(VectorSpace(self.base_ring(), dim))

    def _repr_degrees(self, min_deg, max_deg):
        """
        Return a string representation

        This method is like :meth:`_repr_` except that the user can
        select the range of degrees to be shown in the output.

        INPUT:

        - ``min_deg``, ``max_deg`` -- two integers.

        EXAMPLES::

            sage: rays = [(1,0), (1,1), (1,2), (-1,-1)]
            sage: F = FilteredVectorSpace(rays, {0:[1, 2], 2:[3]})
            sage: F._repr_degrees(-2, 4)
            ['QQ^2', 'QQ^2', 'QQ^2', 'QQ^1', 'QQ^1', '0', '0', '0']
        """
        degrees = list(range(min_deg, max_deg + 1))
        dims = []
        for i in degrees + [infinity]:
            d = self.get_degree(i).dimension()
            dims.append(self._repr_vector_space(d))
        return dims

    def _repr_(self):
        r"""
        Return as string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: rays = [(1,0), (1,1), (1,2), (-1,-1)]
            sage: FilteredVectorSpace(rays, {0:[1, 2], 2:[3]})._repr_()
            'QQ^2 >= QQ^1 >= QQ^1 >= 0'
            sage: FilteredVectorSpace(rays, {0:[1, 2], oo:[3]})
            QQ^2 >= QQ^1
            sage: FilteredVectorSpace(rays, {oo:[3]})
            QQ^1 in QQ^2
            sage: FilteredVectorSpace(rays, {0:[3]})
            QQ^1 >= 0 in QQ^2
            sage: FilteredVectorSpace({1:[(1,0), (-1,1)], 3:[(1,0)]}, base_ring=GF(3))
            GF(3)^2 >= GF(3)^1 >= GF(3)^1 >= 0
            sage: FilteredVectorSpace({1:[(1,0), (-1,1)], 3:[(1,0)]}, base_ring=AA)
            Vector space of dimension 2 over Algebraic Real Field
            >= Vector space of dimension 1 over Algebraic Real Field
            >= Vector space of dimension 1 over Algebraic Real Field >= 0
        """
        finite_support = [d for d in self.support() if d != infinity]
        if len(finite_support) == 0:
            dims = self._repr_degrees(0, -1)
        else:
            min_deg = finite_support[0]
            max_deg = finite_support[-1]
            dims = self._repr_degrees(min_deg, max_deg)
        s = ' >= '.join(dims)
        if not self.is_exhaustive():
            s += ' in ' + self._repr_vector_space(self.degree())
        return s

    def __eq__(self, other):
        """
        Return whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: V = FilteredVectorSpace(2, 0)
            sage: W = FilteredVectorSpace([(1,0),(0,1)], {0:[0, 1]})
            sage: V == W
            True
            sage: V is W
            False

            sage: W = FilteredVectorSpace([(1,0),(1,1)], {0:[1]})
            sage: V == W
            False

        TESTS::

            sage: P = toric_varieties.P2()
            sage: T_P = P.sheaves.tangent_bundle()
            sage: O_P = P.sheaves.trivial_bundle(1)
            sage: S1 = T_P + O_P
            sage: S2 = O_P + T_P
            sage: S1._filt[0].is_isomorphic(S2._filt[0])  # known bug
            True

            sage: FilteredVectorSpace(2, base_ring=QQ) == FilteredVectorSpace(2, base_ring=GF(5))
            False
        """
        if type(self) != type(other):
            return False
        if self.base_ring() != other.base_ring():
            return False
        if self.dimension() != other.dimension():
            return False
        if len(self._filt) != len(other._filt):
            return False
        for self_filt, other_filt in zip(self._filt, other._filt):
            if self_filt[0] != other_filt[0]:
                # compare degree
                return False
            if (self_filt[1].echelonized_basis_matrix() !=
                    other_filt[1].echelonized_basis_matrix()):
                # compare vector subspace
                return False
        return True

    def __ne__(self, other):
        """
        Return whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: V = FilteredVectorSpace(2, 0)
            sage: W = FilteredVectorSpace([(1,0),(0,1)], {0:[0, 1]})
            sage: V != W
            False

            sage: W = FilteredVectorSpace([(1,0),(1,1)], {0:[1]})
            sage: V != W
            True
        """
        return not (self == other)

    def direct_sum(self, other):
        """
        Return the direct sum.

        INPUT:

        - ``other`` -- a filtered vector space.

        OUTPUT:

        The direct sum as a filtered vector space.

        EXAMPLES::

            sage: V = FilteredVectorSpace(2, 0)
            sage: W = FilteredVectorSpace({0:[(1,-1),(2,1)], 1:[(1,1)]})
            sage: V.direct_sum(W)
            QQ^4 >= QQ^1 >= 0
            sage: V + W    # syntactic sugar
            QQ^4 >= QQ^1 >= 0
            sage: V + V == FilteredVectorSpace(4, 0)
            True

            sage: W = FilteredVectorSpace([(1,-1),(2,1)], {1:[0,1], 2:[1]})
            sage: V + W
            QQ^4 >= QQ^2 >= QQ^1 >= 0

        A suitable base ring is chosen if they do not match::

            sage: v = [(1,0), (0,1)]
            sage: F1 = FilteredVectorSpace(v, {0:[0], 1:[1]}, base_ring=QQ)
            sage: F2 = FilteredVectorSpace(v, {0:[0], 1:[1]}, base_ring=RDF)
            sage: F1 + F2
            RDF^4 >= RDF^2 >= 0
        """
        from sage.structure.element import get_coercion_model
        base_ring = get_coercion_model().common_parent(self.base_ring(), other.base_ring())
        # construct the generators
        self_gens, self_filt = self.presentation()
        other_gens, other_filt = other.presentation()
        generators = \
            [ list(v) + [base_ring.zero()]*other.dimension() for v in self_gens  ] + \
            [ [base_ring.zero()]*self.dimension() + list(v)  for v in other_gens ]
        # construct the filtration dictionary
        def join_indices(self_indices, other_indices):
            self_indices = tuple(self_indices)
            other_indices = tuple(i + len(self_gens) for i in other_indices)
            return self_indices + other_indices
        filtration = dict()
        self_indices = set()
        other_indices = set()
        degrees = list(self_filt) + list(other_filt)
        for deg in sorted(set(degrees), reverse=True):
            self_indices.update(self_filt.get(deg, []))
            other_indices.update(other_filt.get(deg, []))
            gens = join_indices(self_indices, other_indices)
            filtration[deg] = gens
        return FilteredVectorSpace(generators, filtration, base_ring=base_ring)

    __add__ = direct_sum

    def tensor_product(self, other):
        r"""
        Return the graded tensor product.

        INPUT:

        - ``other`` -- a filtered vector space.

        OUTPUT:

        The graded tensor product, that is, the tensor product of a
        generator of degree `d_1` with a generator in degree `d_2` has
        degree `d_1 + d_2`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(1, 1)
            sage: F2 = FilteredVectorSpace(1, 2)
            sage: F1.tensor_product(F2)
            QQ^1 >= 0
            sage: F1 * F2
            QQ^1 >= 0

            sage: F1.min_degree()
            1
            sage: F2.min_degree()
            2
            sage: (F1*F2).min_degree()
            3

        A suitable base ring is chosen if they do not match::

            sage: v = [(1,0), (0,1)]
            sage: F1 = FilteredVectorSpace(v, {0:[0], 1:[1]}, base_ring=QQ)
            sage: F2 = FilteredVectorSpace(v, {0:[0], 1:[1]}, base_ring=RDF)
            sage: F1 * F2
            RDF^4 >= RDF^3 >= RDF^1 >= 0
        """
        V = self
        W = other
        from sage.structure.element import get_coercion_model
        base_ring = get_coercion_model().common_parent(V.base_ring(), W.base_ring())
        from sage.modules.tensor_operations import VectorCollection, TensorOperation
        V_generators, V_indices = V.presentation()
        W_generators, W_indices = W.presentation()
        V_coll = VectorCollection(V_generators, base_ring, V.dimension())
        W_coll = VectorCollection(W_generators, base_ring, W.dimension())
        T = TensorOperation([V_coll, W_coll], 'product')

        filtration = dict()
        for V_deg in V.support():
            for W_deg in W.support():
                deg = V_deg + W_deg
                indices = filtration.get(deg, set())
                for i in V_indices[V_deg]:
                    for j in W_indices[W_deg]:
                        i_tensor_j = T.index_map(i, j)
                        indices.add(i_tensor_j)
                filtration[deg] = indices
        return FilteredVectorSpace(T.vectors(), filtration, base_ring=base_ring)

    __mul__ = tensor_product

    def _power_operation(self, n, operation):
        """
        Return tensor power operation.

        INPUT:

        - ``n`` -- integer. the number of factors of ``self``.

        - ``operation`` -- string. See
          :class:`~sage.modules.tensor_operations.TensorOperation` for
          details.

        EXAMPLES::

            sage: F = FilteredVectorSpace(1, 1) + FilteredVectorSpace(1, 2);  F
            QQ^2 >= QQ^1 >= 0
            sage: F._power_operation(2, 'symmetric')
            QQ^3 >= QQ^2 >= QQ^1 >= 0
            sage: F._power_operation(2, 'antisymmetric')
            QQ^1 >= 0
        """
        from sage.modules.tensor_operations import VectorCollection, TensorOperation
        generators, indices = self.presentation()
        V = VectorCollection(generators, self.base_ring(), self.dimension())
        T = TensorOperation([V] * n, operation)

        iters = [self.support()] * n
        filtration = dict()
        from sage.categories.cartesian_product import cartesian_product
        for degrees in cartesian_product(iters):
            deg = sum(degrees)
            filt_deg = filtration.get(deg, set())
            for i in cartesian_product([indices.get(d) for d in degrees]):
                pow_i = T.index_map(*i)
                if pow_i is not None:
                    filt_deg.add(pow_i)
            filtration[deg] = filt_deg
        return FilteredVectorSpace(T.vectors(), filtration, base_ring=self.base_ring())


    def exterior_power(self, n):
        """
        Return the `n`-th graded exterior power.

        INPUT:

        - ``n`` -- integer. Exterior product of how many copies of
          ``self``.

        OUTPUT:

        The graded exterior product, that is, the wedge product of a
        generator of degree `d_1` with a generator in degree `d_2` has
        degree `d_1 + d_2`.

        EXAMPLES::

            sage: F = FilteredVectorSpace(1, 1) + FilteredVectorSpace(1, 2);  F
            QQ^2 >= QQ^1 >= 0
            sage: F.exterior_power(1)
            QQ^2 >= QQ^1 >= 0
            sage: F.exterior_power(2)
            QQ^1 >= 0
            sage: F.exterior_power(3)
            0
            sage: F.wedge(2)
            QQ^1 >= 0
        """
        return self._power_operation(n, 'antisymmetric')

    wedge = exterior_power

    def symmetric_power(self, n):
        """
        Return the `n`-th graded symmetric power.

        INPUT:

        - ``n`` -- integer. Symmetric product of how many copies of
          ``self``.

        OUTPUT:

        The graded symmetric product, that is, the symmetrization of a
        generator of degree `d_1` with a generator in degree `d_2` has
        degree `d_1 + d_2`.

        EXAMPLES::

            sage: F = FilteredVectorSpace(1, 1) + FilteredVectorSpace(1, 2);  F
            QQ^2 >= QQ^1 >= 0
            sage: F.symmetric_power(2)
            QQ^3 >= QQ^2 >= QQ^1 >= 0
        """
        return self._power_operation(n, 'symmetric')

    def dual(self):
        """
        Return the dual filtered vector space.

        OUTPUT:

        The graded dual, that is, the dual of a degree-`d` subspace is
        a set of linear constraints in degree `-d+1`. That is, the
        dual generators live in degree `-d`.

        EXAMPLES::

            sage: gens = identity_matrix(3).rows()
            sage: F = FilteredVectorSpace(gens, {0:[0,1,2], 2:[0]});  F
            QQ^3 >= QQ^1 >= QQ^1 >= 0
            sage: F.support()
            (0, 2)

            sage: F.dual()
            QQ^3 >= QQ^2 >= QQ^2 >= 0
            sage: F.dual().support()
            (-2, 0)
        """
        filtration = dict()
        prev_deg = minus_infinity
        for deg, V in self._filt[1:]:
            filtration[-prev_deg] = V.complement().echelonized_basis()
            prev_deg = deg
        return FilteredVectorSpace(filtration, base_ring=self.base_ring())

    def shift(self, deg):
        """
        Return a filtered vector space with degrees shifted by a constant.

        EXAMPLES::

            sage: gens = identity_matrix(3).rows()
            sage: F = FilteredVectorSpace(gens, {0:[0,1,2], 2:[0]});  F
            QQ^3 >= QQ^1 >= QQ^1 >= 0
            sage: F.support()
            (0, 2)
            sage: F.shift(-5).support()
            (-5, -3)
        """
        generators, filtration = self.presentation()
        shifted = dict()
        for d, indices in filtration.items():
            shifted[d + deg] = indices
        return FilteredVectorSpace(generators, shifted, base_ring=self.base_ring())

    def random_deformation(self, epsilon=None):
        """
        Return a random deformation

        INPUT:

        - ``epsilon`` -- a number in the base ring.

        OUTPUT:

        A new filtered vector space where the generators of the
        subspaces are moved by ``epsilon`` times a random vector.

        EXAMPLES::

            sage: gens = identity_matrix(3).rows()
            sage: F = FilteredVectorSpace(gens, {0:[0,1,2], 2:[0]});  F
            QQ^3 >= QQ^1 >= QQ^1 >= 0
            sage: F.get_degree(2)
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [1 0 0]
            sage: G = F.random_deformation(1/50);  G
            QQ^3 >= QQ^1 >= QQ^1 >= 0
            sage: D = G.get_degree(2)
            sage: D.degree()
            3
            sage: v = D.basis_matrix()[0]
            sage: v[0]
            1

            sage: while F.random_deformation(1/50).get_degree(2).matrix() == matrix([1, 0, 0]):
            ....:     pass
        """
        from sage.modules.free_module_element import random_vector
        R = self.base_ring()
        if epsilon is None:
            epsilon = R.one()
        filtration = dict()
        for deg, filt in self._filt[1:]:
            generators = [v + epsilon * random_vector(R, self.rank())
                          for v in filt.echelonized_basis()]
            filtration[deg] = generators
        return FilteredVectorSpace(filtration, base_ring=R, check=True)
