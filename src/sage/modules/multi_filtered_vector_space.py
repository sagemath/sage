r"""
Multiple `\ZZ`-Graded Filtrations of a Single Vector Space

See :mod:`filtered_vector_space` for simply graded vector spaces. This
module implements the analog but for a collection of filtrations of
the same vector space.

The basic syntax to use it is a dictionary whose keys are some
arbitrary indexing set and values are
:func:`~sage.modules.filtered_vector_space.FilteredVectorSpace` ::

    sage: F1 = FilteredVectorSpace(2, 1)
    sage: F2 = FilteredVectorSpace({0:[(1,0)], 2:[(2,3)]})
    sage: V = MultiFilteredVectorSpace({'first':F1, 'second':F2})
    sage: V
    Filtrations
         first: QQ^2 >= QQ^2 >=  0   >= 0
        second: QQ^2 >= QQ^1 >= QQ^1 >= 0

    sage: V.index_set()   # random output
    {'second', 'first'}
    sage: sorted(V.index_set())
    ['first', 'second']

    sage: V.get_filtration('first')
    QQ^2 >=  0
    sage: V.get_degree('second', 1)
    Vector space of degree 2 and dimension 1 over Rational Field
    Basis matrix:
    [  1 3/2]
"""

# ****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.all import QQ, ZZ, Integer
from sage.rings.infinity import infinity, minus_infinity
from sage.categories.fields import Fields
from sage.modules.free_module import FreeModule_ambient_field, VectorSpace
from sage.misc.cachefunc import cached_method
from sage.modules.filtered_vector_space import FilteredVectorSpace


def MultiFilteredVectorSpace(arg, base_ring=None, check=True):
    """
    Contstruct a multi-filtered vector space.

    INPUT:

    - ``arg`` -- either a non-empty dictionary of filtrations or an
      integer. The latter is interpreted as the vector space
      dimension, and the indexing set of the filtrations is empty.

    - ``base_ring`` -- a field (optional, default ``'None'``). The
      base field of the vector space. Must be a field. If not
      specified, the base field is derived from the filtrations.

    - ``check`` -- boolean (optional; default: ``True``). Whether
      to perform consistency checks.

    EXAMPLES::

        sage: MultiFilteredVectorSpace(3, QQ)
        Unfiltered QQ^3

        sage: F1 = FilteredVectorSpace(2, 1)
        sage: F2 = FilteredVectorSpace(2, 3)
        sage: V = MultiFilteredVectorSpace({1:F1, 2:F2});  V
        Filtrations
            1: QQ^2 >=  0   >=  0   >= 0
            2: QQ^2 >= QQ^2 >= QQ^2 >= 0
   """
    if arg in ZZ:
        dim = ZZ(arg)
        filtration = {}
        if base_ring is None:
            base_ring = QQ
    else:
        filtration = dict(arg)
        F = next(iter(arg.values()))  # the first filtration
        dim = F.dimension()
        if base_ring is None:
            base_ring = F.base_ring()
    for deg in filtration:
        filt = filtration[deg]
        if filt.base_ring() != base_ring:
            filt = filt.change_ring(base_ring)
        filtration[deg] = filt
    return MultiFilteredVectorSpace_class(base_ring, dim, filtration)


class MultiFilteredVectorSpace_class(FreeModule_ambient_field):

    def __init__(self, base_ring, dim, filtrations, check=True):
        """
        Python constructor.

        .. warning::

            Use :func:`MultiFilteredVectorSpace` to construct
            multi-filtered vector spaces.

        INPUT:

        - ``base_ring`` -- a ring. the base ring.

        - ``dim`` -- integer. The dimension of the ambient vector space.

        - ``filtrations`` -- a dictionary whose values are
          filtrations.

        - ``check`` -- boolean (optional). Whether to perform
          additional consistency checks.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2});  V
            Filtrations
                1: QQ^2 >=  0   >=  0   >= 0
                2: QQ^2 >= QQ^2 >= QQ^2 >= 0
        """
        if check:
            assert isinstance(dim, Integer)
            assert base_ring in Fields()
            assert all(base_ring == f.base_ring() for f in filtrations.values())
            assert all(dim == f.dimension() for f in filtrations.values())
        super(MultiFilteredVectorSpace_class, self).__init__(base_ring, dim)
        self._filt = dict(filtrations)

    @cached_method
    def index_set(self):
        """
        Return the allowed indices for the different filtrations.

        OUTPUT:

        Set.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.index_set()
            {1, 2}
        """
        from sage.sets.set import Set
        return Set(self._filt.keys())

    def change_ring(self, base_ring):
        """
        Return the same multi-filtration over a different base ring.

        INPUT:

        - ``base_ring`` -- a ring. The new base ring.

        OUTPUT:

        This method returns a new multi-filtered vector space whose
        subspaces are defined by the same generators but over a
        different base ring.

        EXAMPLES::

            sage: V = FilteredVectorSpace(2, 0)
            sage: W = FilteredVectorSpace(2, 2)
            sage: F = MultiFilteredVectorSpace({'a':V, 'b':W});  F
            Filtrations
                a: QQ^2 >=  0   >=  0   >= 0
                b: QQ^2 >= QQ^2 >= QQ^2 >= 0
            sage: F.change_ring(RDF)
            Filtrations
                a: RDF^2 >=   0   >=   0   >= 0
                b: RDF^2 >= RDF^2 >= RDF^2 >= 0

            sage: MultiFilteredVectorSpace(3, base_ring=QQ).change_ring(RR)
            Unfiltered RR^3
        """
        if not self._filt:
            return MultiFilteredVectorSpace(self.dimension(),
                                            base_ring=base_ring)
        filtrations = {}
        for key, F in self._filt.items():
            filtrations[key] = F.change_ring(base_ring)
        return MultiFilteredVectorSpace(filtrations, base_ring=base_ring)

    def ambient_vector_space(self):
        """
        Return the ambient (unfiltered) vector space.

        OUTPUT:

        A vector space.

        EXAMPLES::

            sage: V = FilteredVectorSpace(2, 0)
            sage: W = FilteredVectorSpace(2, 2)
            sage: F = MultiFilteredVectorSpace({'a':V, 'b':W})
            sage: F.ambient_vector_space()
            Vector space of dimension 2 over Rational Field
        """
        return VectorSpace(self.base_ring(), self.dimension())

    @cached_method
    def is_constant(self):
        """
        Return whether the multi-filtration is constant.

        OUTPUT:

        Boolean. Whether the each filtration is constant, see
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.is_constant`.

        EXAMPLES::

            sage: V = FilteredVectorSpace(2, 0)
            sage: W = FilteredVectorSpace(2, 2)
            sage: F = MultiFilteredVectorSpace({'a':V, 'b':W});  F
            Filtrations
                a: QQ^2 >=  0   >=  0   >= 0
                b: QQ^2 >= QQ^2 >= QQ^2 >= 0
            sage: F.is_constant()
            False
        """
        return all(F.is_constant() for F in self._filt.values())

    def is_exhaustive(self):
        r"""
        Return whether the multi-filtration is exhaustive.

        A filtration $\{F_d\}$ in an ambient vector space $V$ is
        exhaustive if $\cup F_d = V$. See also :meth:`is_separating`.

        OUTPUT:

        Boolean. Whether each filtration is constant, see
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.is_exhaustive`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.is_exhaustive()
            True
        """
        return all(F.is_exhaustive() for F in self._filt.values())

    def is_separating(self):
        r"""
        Return whether the multi-filtration is separating.

        A filtration $\{F_d\}$ in an ambient vector space $V$ is
        exhaustive if $\cap F_d = 0$. See also :meth:`is_exhaustive`.

        OUTPUT:

        Boolean. Whether each filtration is separating, see
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.is_separating`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.is_separating()
            True
        """
        return all(F.is_separating() for F in self._filt.values())

    @cached_method
    def support(self):
        """
        Return the degrees in which there are non-trivial generators.

        OUTPUT:

        A tuple of integers (and plus infinity) in ascending
        order. The last entry is plus infinity if and only if the
        filtration is not separating (see :meth:`is_separating`).

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.support()
            (1, 3)
        """
        support = set()
        for F in self._filt.values():
            support.update(F.support())
        return tuple(sorted(support))

    @cached_method
    def min_degree(self):
        r"""
        Return the lowest degree of the filtration.

        OUTPUT:

        Integer or plus infinity. The largest degree `d` of the
        (descending) filtrations such that, for each individual
        filtration, the filtered vector space `F_d` still equal to
        `F_{-\infty}`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.min_degree()
            1
        """
        if not self._filt:
            return infinity
        return min(F.min_degree() for F in self._filt.values())

    @cached_method
    def max_degree(self):
        r"""
        Return the highest degree of the filtration.

        OUTPUT:

        Integer or minus infinity. The smallest degree of the
        filtrations such that the filtrations are constant to the
        right.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.max_degree()
            4
        """
        if not self._filt:
            return minus_infinity
        return max(F.max_degree() for F in self._filt.values())

    def get_filtration(self, key):
        """
        Return the filtration indexed by ``key``.

        OUTPUT:

        A filtered vector space.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.get_filtration(2)
            QQ^2 >= 0
        """
        return self._filt[key]

    def get_degree(self, key, deg):
        r"""
        Return one filtered vector space.

        INPUT:

        - ``key`` -- an element of the :meth:`index_set`. Specifies
          which filtration.

        - ``d`` -- Integer. The desired degree of the filtration.

        OUTPUT:

        The vector space of degree ``deg`` in the filtration indexed
        by ``key`` as subspace of the ambient space.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(2, 3)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.get_degree(2, 0)
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
        """
        return self._filt[key].get_degree(deg)

    def graded(self, key, deg):
        r"""
        Return the associated graded vector space.

        INPUT:

        - ``key`` -- an element of the :meth:`index_set`. Specifies
          which filtration.

        - ``d`` -- Integer. The desired degree of the filtration.

        OUTPUT:

        The quotient `G_d = F_d / F_{d+1}` of the filtration `F`
        corresponding to ``key``.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V.graded(2, 3)
            Vector space quotient V/W of dimension 1 over Rational Field where
            V: Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
            W: Vector space of degree 2 and dimension 0 over Rational Field
            Basis matrix:
            []
        """
        return self.get_degree(key, deg).quotient(self.get_degree(key, deg + 1))

    def _repr_(self):
        r"""
        Return as string representation of ``self``.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: rays = [(1,0), (1,1), (1,2), (-1,-1)]
            sage: F1 = FilteredVectorSpace(rays, {0:[1, 2], 2:[3]})
            sage: F2 = FilteredVectorSpace(rays, {0:[1, 2], oo:[3]})
            sage: F3 = FilteredVectorSpace(rays, {oo:[3]})
            sage: F4 = FilteredVectorSpace(rays, {0:[3]})
            sage: MultiFilteredVectorSpace({'a':F1, 'b':F2, 'c': F3, 'd': F4})
            Filtrations
                a: QQ^2 >= QQ^1 >= QQ^1 >=  0
                b: QQ^2 >= QQ^1 >= QQ^1 >= QQ^1
                c: QQ^1 >= QQ^1 >= QQ^1 >= QQ^1
                d: QQ^1 >=  0   >=  0   >=  0

            sage: MultiFilteredVectorSpace(123, base_ring=RR)
            Unfiltered RR^123
        """
        if not self._filt:
            F = FilteredVectorSpace(self.dimension(),
                                    base_ring=self.base_ring())
            return 'Unfiltered ' + repr(F)
        rows = []
        min_deg, max_deg = self.min_degree(), self.max_degree()
        for key in sorted(self.index_set()):
            F = self.get_filtration(key)
            r = [str(key)] + F._repr_degrees(min_deg, max_deg - 1)
            rows.append(r)
        from sage.misc.table import table
        t = table(rows)
        w = t._widths()
        lines = ['Filtrations']
        for r in rows:
            s = '    '
            s += r[0].rjust(w[0]) + ': '
            s += ' >= '.join(r[i].center(w[i]) for i in range(1, len(w)))
            lines.append(s)
        return '\n'.join(lines)

    def __eq__(self, other):
        """
        Return whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V == MultiFilteredVectorSpace({2:F2, 1:F1})
            True
            sage: V == MultiFilteredVectorSpace({'a':F1, 'b':F2})
            False
        """
        if type(self) != type(other):
            return False
        return self._filt == other._filt

    def __ne__(self, other):
        """
        Return whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({1:F1, 2:F2})
            sage: V != MultiFilteredVectorSpace({2:F2, 1:F1})
            False
            sage: V != MultiFilteredVectorSpace({'a':F1, 'b':F2})
            True
        """
        return not (self == other)

    def direct_sum(self, other):
        """
        Return the direct sum.

        INPUT:

        - ``other`` -- a multi-filtered vector space with the same
          :meth:`index_set`.

        OUTPUT:

        The direct sum as a multi-filtered vector space. See
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.direct_sum`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({'a':F1, 'b':F2})
            sage: G1 = FilteredVectorSpace(1, 1)
            sage: G2 = FilteredVectorSpace(1, 3)
            sage: W = MultiFilteredVectorSpace({'a':G1, 'b':G2})
            sage: V.direct_sum(W)
            Filtrations
                a: QQ^3 >= QQ^3 >=  0   >=  0   >= 0
                b: QQ^3 >= QQ^2 >= QQ^2 >= QQ^2 >= 0
            sage: V + W   # syntactic sugar
            Filtrations
                a: QQ^3 >= QQ^3 >=  0   >=  0   >= 0
                b: QQ^3 >= QQ^2 >= QQ^2 >= QQ^2 >= 0
        """
        if not self.index_set() == other.index_set():
            raise ValueError('the index sets of the two summands'
                             ' must be the same')
        filtrations = {}
        for key in self.index_set():
            filtrations[key] = self._filt[key] + other._filt[key]
        return MultiFilteredVectorSpace(filtrations)

    __add__ = direct_sum

    def tensor_product(self, other):
        r"""
        Return the graded tensor product.

        INPUT:

        - ``other`` -- a multi-filtered vector space with the same
          :meth:`index_set`.

        OUTPUT:

        The tensor product of ``self`` and ``other`` as a
        multi-filtered vector space. See
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.tensor_product`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({'a':F1, 'b':F2})
            sage: G1 = FilteredVectorSpace(1, 1)
            sage: G2 = FilteredVectorSpace(1, 3)
            sage: W = MultiFilteredVectorSpace({'a':G1, 'b':G2})
            sage: V.tensor_product(W)
            Filtrations
                a: QQ^2 >=  0   >=  0   >=  0   >=  0   >= 0
                b: QQ^2 >= QQ^2 >= QQ^1 >= QQ^1 >= QQ^1 >= 0
            sage: V * W   # syntactic sugar
            Filtrations
                a: QQ^2 >=  0   >=  0   >=  0   >=  0   >= 0
                b: QQ^2 >= QQ^2 >= QQ^1 >= QQ^1 >= QQ^1 >= 0
        """
        if not self.index_set() == other.index_set():
            raise ValueError('the index sets of the two summands'
                             ' must be the same')
        filtrations = {}
        for key in self.index_set():
            filtrations[key] = self._filt[key] * other._filt[key]
        return MultiFilteredVectorSpace(filtrations)

    __mul__ = tensor_product

    def exterior_power(self, n):
        """
        Return the `n`-th graded exterior power.

        INPUT:

        - ``n`` -- integer. Exterior product of how many copies of
          ``self``.

        OUTPUT:

        The exterior power as a multi-filtered vector space. See
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.exterior_power`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({'a':F1, 'b':F2})
            sage: V.exterior_power(2)  # long time
            Filtrations
                a: QQ^1 >=  0   >= 0
                b: QQ^1 >= QQ^1 >= 0
        """
        filtrations = {key: value.exterior_power(n)
                       for key, value in self._filt.items()}
        return MultiFilteredVectorSpace(filtrations)

    wedge = exterior_power

    def symmetric_power(self, n):
        """
        Return the `n`-th graded symmetric power.

        INPUT:

        - ``n`` -- integer. Symmetric product of how many copies of
          ``self``.

        OUTPUT:

        The symmetric power as a multi-filtered vector space. See
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.symmetric_power`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({'a':F1, 'b':F2})
            sage: V.symmetric_power(2)
            Filtrations
                a: QQ^3 >= QQ^3 >= QQ^3 >=  0   >=  0   >=  0   >=  0   >= 0
                b: QQ^3 >= QQ^2 >= QQ^2 >= QQ^2 >= QQ^1 >= QQ^1 >= QQ^1 >= 0
        """
        filtrations = {key: value.symmetric_power(n)
                       for key, value in self._filt.items()}
        return MultiFilteredVectorSpace(filtrations)

    def dual(self):
        """
        Return the dual.

        OUTPUT:

        The dual as a multi-filtered vector space. See
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.dual`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({'a':F1, 'b':F2})
            sage: V.dual()
            Filtrations
                a: QQ^2 >= QQ^2 >= QQ^2 >=  0   >= 0
                b: QQ^2 >= QQ^1 >= QQ^1 >= QQ^1 >= 0
        """
        filtrations = {key: value.dual()
                       for key, value in self._filt.items()}
        return MultiFilteredVectorSpace(filtrations)

    def shift(self, deg):
        """
        Return a filtered vector space with degrees shifted by a constant.

        OUTPUT:

        The shift of ``self``. See
        :meth:`~sage.modules.filtered_vector_space.FilteredVectorSpace_class.shift`.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({'a':F1, 'b':F2})
            sage: V.support()
            (0, 1, 3)
            sage: V.shift(-5).support()
            (-5, -4, -2)
        """
        filtrations = {key: value.shift(deg)
                       for key, value in self._filt.items()}
        return MultiFilteredVectorSpace(filtrations)

    def random_deformation(self, epsilon=None):
        """
        Return a random deformation

        INPUT:

        - ``epsilon`` -- a number in the base ring.

        OUTPUT:

        A new multi-filtered vector space where the generating vectors
        of subspaces are moved by ``epsilon`` times a random vector.

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(2, 1)
            sage: F2 = FilteredVectorSpace(1, 3) + FilteredVectorSpace(1,0)
            sage: V = MultiFilteredVectorSpace({'a':F1, 'b':F2})
            sage: V.get_degree('b',1)
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [1 0]
            sage: D = V.random_deformation(1/100).get_degree('b',1)
            sage: D.degree()
            2
            sage: D.dimension()
            1
            sage: D.matrix()[0, 0]
            1

            sage: while V.random_deformation(1/100).get_degree('b',1).matrix() == matrix([1, 0]):
            ....:     pass
        """
        filtrations = {key: value.random_deformation(epsilon)
                       for key, value in self._filt.items()}
        return MultiFilteredVectorSpace(filtrations)
