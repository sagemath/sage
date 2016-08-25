r"""
Point collections

This module was designed as a part of framework for toric varieties
(:mod:`~sage.schemes.toric.variety`,
:mod:`~sage.schemes.toric.fano_variety`).

AUTHORS:

- Andrey Novoseltsev (2011-04-25): initial version, based on cone module.

- Andrey Novoseltsev (2012-03-06): additions and doctest changes while
  switching cones to use point collections.

EXAMPLES:

The idea behind :class:`point collections <PointCollection>` is to have a
container for points of the same space that

* behaves like a tuple *without significant performance penalty*::

    sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
    sage: c[1]
    N(1, 0, 1)
    sage: for point in c: point
    N(0, 0, 1)
    N(1, 0, 1)
    N(0, 1, 1)
    N(1, 1, 1)

* prints in a convenient way and with clear indication of the ambient space::

    sage: c
    N(0, 0, 1),
    N(1, 0, 1),
    N(0, 1, 1),
    N(1, 1, 1)
    in 3-d lattice N

* allows (cached) access to alternative representations::

    sage: c.set()
    frozenset({N(0, 0, 1), N(0, 1, 1), N(1, 0, 1), N(1, 1, 1)})

* allows introduction of additional methods::

    sage: c.basis()
    N(0, 0, 1),
    N(1, 0, 1),
    N(0, 1, 1)
    in 3-d lattice N

Examples of natural point collections include ray and line generators of cones,
vertices and points of polytopes, normals to facets, their subcollections, etc.

Using this class for all of the above cases allows for unified interface *and*
cache sharing. Suppose that `\Delta` is a reflexive polytope. Then the same
point collection can be linked as

#. vertices of `\Delta`;
#. facet normals of its polar `\Delta^\circ`;
#. ray generators of the face fan of `\Delta`;
#. ray generators of the normal fan of `\Delta`.

If all these objects are in use and, say, a matrix representation was computed
for one of them, it becomes available to all others as well, eliminating the
need to spend time and memory four times.
"""

#*****************************************************************************
#       Copyright (C) 2012 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.structure.sage_object cimport SageObject

from sage.matrix.all import matrix
from sage.misc.all import latex


def is_PointCollection(x):
    r"""
    Check if ``x`` is a :class:`point collection <PointCollection>`.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a point collection and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.point_collection import is_PointCollection
        sage: is_PointCollection(1)
        False
        sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
        sage: is_PointCollection(c.rays())
        True
    """
    return isinstance(x, PointCollection)


_output_format = "default"


cdef class PointCollection(SageObject):
    r"""
    Create a point collection.

    .. WARNING::

        No correctness check or normalization is performed on the input data.
        This class is designed for internal operations and you probably should
        not use it directly.

    Point collections are immutable, but cache most of the returned values.

    INPUT:

    - ``points`` -- an iterable structure of immutable elements of ``module``,
      if ``points`` are already accessible to you as a :class:`tuple`, it is
      preferable to use it for speed and memory consumption reasons;

    - ``module`` -- an ambient module for ``points``. If ``None``, it will be
      determined as :func:`parent` of the first point. Of course, this cannot
      be done if there are no points, so in this case you must give an
      appropriate ``module`` directly. Note that ``None`` is *not* the default
      value - you always *must* give this argument explicitly, even if it is
      ``None``.

    OUTPUT:

    - a point collection.
    """

    cdef:
        tuple _points
        object _module
        # cache attributes
        PointCollection _basis
        object _matrix
        frozenset _set

    def __init__(self, points, module=None):
        r"""
        See :class:`PointCollection` for documentation.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: v = vector([1,0])
            sage: v.set_immutable()
            sage: c = PointCollection([v], ZZ^2)
            sage: c.module()
            Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: c = PointCollection([v], None)
            sage: c.module()  # Determined automatically
            Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: TestSuite(c).run()
        """
        super(PointCollection, self).__init__()
        self._points = tuple(points)
        self._module = self._points[0].parent() if module is None else module

    def __add__(left, right):
        r"""
        Return the joint point collection.

        INPUT:

        - ``left`` -- a :class:`PointCollection`;

        - ``right`` -- a :class:`PointCollection`.

        OUTPUT:

        - a :class:`PointCollection`.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c + c
            N(0, 0, 1),
            N(1, 0, 1),
            N(0, 1, 1),
            N(1, 1, 1),
            N(0, 0, 1),
            N(1, 0, 1),
            N(0, 1, 1),
            N(1, 1, 1)
            in 3-d lattice N
        """
        if not (isinstance(left, PointCollection) and
                isinstance(right, PointCollection)):
            raise NotImplementedError
        cdef PointCollection left_pc = left
        cdef PointCollection right_pc = right
        if not left_pc._module is right_pc._module:
            raise NotImplementedError
        return PointCollection(left_pc._points + right_pc._points,
                               left_pc._module)

    def __call__(self, *args):
        r"""
        Return a subcollection of ``self``.

        INPUT:

        - a list of integers (as a single or many arguments).

        OUTPUT:

        - a :class:`point collection <PointCollection>`.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c()
            Empty collection
            in 3-d lattice N
            sage: c(2,1)
            N(0, 1, 1),
            N(1, 0, 1)
            in 3-d lattice N
            sage: c(range(4)) == c
            True
        """
        if len(args) == 1:
            try:
                args = tuple(args[0])
            except TypeError:
                pass
        # Avoid creating a copy of self
        if len(args) == len(self) and args == tuple(range(len(self))):
            return self
        else:
            return PointCollection([self[i] for i in args], self._module)

    def __cmp__(self, right):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        - 0 if ``right`` is of the same type as ``self`` (i.e. it is another
          :class:`point collection <PointCollection>`), they have the same
          :meth:`module`, and their points are the same and listed in the same
          order. 1 or -1 otherwise.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: cmp(c, c)
            0
            sage: cmp(c, 1) * cmp(1, c)
            -1
        """
        c = cmp(type(self), type(right))
        if c != 0:
            return c
        cdef PointCollection pc = right
        c = cmp(self._module, pc._module)
        if c != 0:
            return c
        return cmp(self._points, pc._points)

    def __getitem__(self, n):
        r"""
        Return the ``n``-th point of ``self``.

        INPUT:

        - ``n`` -- an integer.

        OUTPUT:

        - a point, an element of the ambient :meth:`module` of ``self``.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c[0]
            N(0, 0, 1)
        """
        return self._points[n]

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        OUTPUT:

        - an integer.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: hash(c) == hash(c)
            True
        """
        return hash(self._points)

    def __iter__(self):
        r"""
        Return an iterator over points of ``self``.

        OUTPUT:

        - an iterator.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: for point in c: print(point)
            N(0, 0, 1)
            N(1, 0, 1)
            N(0, 1, 1)
            N(1, 1, 1)
        """
        return iter(self._points)

    def __len__(self):
        r"""
        Return the number of points in ``self``.

        OUTPUT:

        - an integer.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: len(c)
            4
        """
        return len(self._points)

    def __list__(self):
        r"""
        Return a list of points of ``self``.

        OUTPUT:

        - a list.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: list(c)
            [N(0, 0, 1), N(1, 0, 1), N(0, 1, 1), N(1, 1, 1)]
        """
        return list(self._points)

    def __mul__(left, right):
        r"""
        Return the product ``left * right``.

        INPUT:

        - a :class:`point collection <PointCollection>` and something that can
          act both on ``self.module().zero()`` and either ``self.matrix()`` from
          the right or ``self.column_matrix()`` from the left.

        OUTPUT:

        - the result of ``self.matrix() * right``, provided that
          ``self.module().zero() * right`` can be computed.

        The idea of this method is to provide a shortcut for matrix
        multiplication with appropriate type checks, in particular, it is not
        possible to multiply by a point of the same toric lattice as elements of
        ``self``.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.matrix()
            [0 0 1]
            [1 0 1]
            [0 1 1]
            [1 1 1]
            sage: c * c[0]
            Traceback (most recent call last):
            ...
            TypeError: elements of the same toric lattice cannot be multiplied!

        If you really need such a product, state it explicitly::

            sage: c.matrix() * c[0]
            (1, 1, 1, 1)

        Multiplication by matrices works as well::

            sage: c * c.column_matrix()
            [1 1 1 1]
            [1 2 1 2]
            [1 1 2 2]
            [1 2 2 3]
        """
        cdef PointCollection pc
        if isinstance(left, PointCollection):
            pc = left
            # Check that it is possible to act on points.
            pc._module.zero() * right
            return pc.matrix() * right
        if isinstance(right, PointCollection):
            pc = right
            # Check that it is possible to act on points.
            left * pc._module.zero()
            return left * pc.column_matrix()
        raise NotImplementedError

    def __reduce__(self):
        r"""
        Prepare ``self`` for pickling.

        OUTPUT:

        - a tuple, currently the class name and a tuple consisting of points
          and the ambient module.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: loads(dumps(c))
            N(0, 0, 1),
            N(1, 0, 1),
            N(0, 1, 1),
            N(1, 1, 1)
            in 3-d lattice N
            sage: loads(dumps(c)) == c
            True
        """
        return (PointCollection, (self._points, self._module))

    def __tuple__(self):
        r"""
        Return the tuple of points of ``self``.

        OUTPUT:

        - a tuple.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: tuple(c)
            (N(0, 0, 1), N(1, 0, 1), N(0, 1, 1), N(1, 1, 1))
        """
        return self._points

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: print(c._latex_())
            \left(\left(0,\,0,\,1\right)_{N}, \left(1,\,0,\,1\right)_{N},
            \left(0,\,1,\,1\right)_{N}, \left(1,\,1,\,1\right)_{N}\right)_{N}
        """
        global _output_format
        if _output_format in ["default", "tuple"]:
            r = latex(tuple(self))
        elif _output_format == "matrix":
            r = latex(self.matrix())
        elif _output_format == "column matrix":
            r = latex(self.column_matrix())
        elif _output_format == "separated column matrix":
            r = latex(self.column_matrix())
            r = r.replace("r" * len(self), "|".join("r" * len(self)))
        return r"%s_{%s}" % (r, latex(self.module()))

    def _matrix_(self, ring=None):
        r"""
        Return a matrix whose rows are points of ``self``.

        INPUT:

        - ``ring`` -- a base ring for the returned matrix (default: base ring of
          :meth:`module` of ``self``).

        OUTPUT:

        - a :class:`matrix <Matrix>`.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: matrix(c) # indirect doctest
            [0 0 1]
            [1 0 1]
            [0 1 1]
            [1 1 1]
        """
        if ring is None:
            return self.matrix()
        else:
            return self.matrix().change_ring(ring)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: print(c._repr_())
            N(0, 0, 1),
            N(1, 0, 1),
            N(0, 1, 1),
            N(1, 1, 1)
            in 3-d lattice N
        """
        global _output_format
        if _output_format == "default":
            r = map(repr, self)
            r = [point.split(",") for point in r]
            if not r:
                r = "Empty collection"
            else:
                if "(" in r[0][0]:
                    delimiter = "("
                elif "[" in r[0][0]:
                    delimiter = "["
                else:
                    raise ValueError("cannot parse point representation!")
                heads = []
                for point in r:
                    head, point[0] = point[0].rsplit(delimiter, 1)
                    heads.append(head + delimiter)
                format = "{{:<{}}}".format(max(map(len, heads)))
                widths = [0] * len(r[0])
                for point in r:
                    for i, coordinate in enumerate(point):
                        widths[i] = max(widths[i], len(coordinate))
                format += ",".join("{{:>{}}}".format(width) for width in widths)
                r = ",\n".join([format.format(head, *point)
                                for head, point in zip(heads, r)])
        elif _output_format == "tuple":
            r = tuple(self)
        elif _output_format == "matrix":
            r = self.matrix()
        else:
            r = self.column_matrix()
        return "{}\nin {}".format(r, self.module())

    def basis(self):
        r"""
        Return a linearly independent subset of points of ``self``.

        OUTPUT:

        - a :class:`point collection <PointCollection>` giving a random (but
          fixed) choice of an `\RR`-basis for the vector space spanned by the
          points of ``self``.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.basis()
            N(0, 0, 1),
            N(1, 0, 1),
            N(0, 1, 1)
            in 3-d lattice N

        Calling this method twice will always return *exactly the same* point
        collection::

            sage: c.basis().basis() is c.basis()
            True
        """
        if self._basis is None:
            self._basis = self(self.matrix().pivot_rows())
        return self._basis

    def cardinality(self):
        r"""
        Return the number of points in ``self``.

        OUTPUT:

        - an integer.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.cardinality()
            4
        """
        return len(self._points)

    def cartesian_product(self, other, module=None):
        r"""
        Return the Cartesian product of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a :class:`point collection <PointCollection>`;

        - ``module`` -- (optional) the ambient module for the result. By
          default, the direct sum of the ambient modules of ``self`` and
          ``other`` is constructed.

        OUTPUT:

        - a :class:`point collection <PointCollection>`.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,1,1)]).rays()
            sage: c.cartesian_product(c)
            N+N(0, 0, 1, 0, 0, 1),
            N+N(1, 1, 1, 0, 0, 1),
            N+N(0, 0, 1, 1, 1, 1),
            N+N(1, 1, 1, 1, 1, 1)
            in 6-d lattice N+N
        """
        assert is_PointCollection(other)
        if module is None:
            module = self._module.direct_sum(other.module())
        P = [list(p) for p in self]
        Q = [list(q) for q in other]
        PQ = [module(p + q) for q in Q for p in P]
        for pq in PQ:
            pq.set_immutable()
        return PointCollection(PQ, module)

    def column_matrix(self):
        r"""
        Return a matrix whose columns are points of ``self``.

        OUTPUT:

        - a :class:`matrix <Matrix>`.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.column_matrix()
            [0 1 0 1]
            [0 0 1 1]
            [1 1 1 1]
        """
        return self.matrix().transpose()

    def dimension(self):
        r"""
        Return the dimension of the space spanned by points of ``self``.

        .. note:: You can use either :meth:`dim` or :meth:`dimension`.

        OUTPUT:

        - an integer.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,1,1)]).rays()
            sage: c.dimension()
            2
            sage: c.dim()
            2
        """
        return self.matrix().rank()

    dim = dimension

    def dual_module(self):
        r"""
        Return the dual of the ambient module of ``self``.

        OUTPUT:

        - a :class:`module <FreeModule_generic>`. If possible (that is, if the
          ambient :meth:`module` `M` of ``self`` has a ``dual()`` method), the
          dual module is returned. Otherwise, `R^n` is returned, where `n` is
          the dimension of `M` and `R` is its base ring.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.dual_module()
            3-d lattice M
        """
        M = self._module
        try:
            return M.dual()
        except AttributeError:
            # TODO: add support for torsion modules as well?
            return M.base_ring() ** M.dimension()

    def index(self, *args):
        r"""
        Return the index of the first occurrence of ``point`` in ``self``.

        INPUT:

        - ``point`` -- a point of ``self``;

        - ``start`` -- (optional) an integer, if given, the search will start
          at this position;

        - ``stop`` -- (optional) an integer, if given, the search will stop
          at this position.

        OUTPUT:

        - an integer if ``point`` is in ``self[start:stop]``, otherwise a
          ``ValueError`` exception is raised.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.index((0,1,1))
            Traceback (most recent call last):
            ...
            ValueError: tuple.index(x): x not in tuple

        Note that this was not a mistake: the *tuple* ``(0,1,1)`` is *not* a
        point of ``c``! We need to pass actual element of the ambient module of
        ``c`` to get their indices::

            sage: N = c.module()
            sage: c.index(N(0,1,1))
            2
            sage: c[2]
            N(0, 1, 1)
        """
        return self._points.index(*args)

    def matrix(self):
        r"""
        Return a matrix whose rows are points of ``self``.

        OUTPUT:

        - a :class:`matrix <Matrix>`.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.matrix()
            [0 0 1]
            [1 0 1]
            [0 1 1]
            [1 1 1]
        """
        if self._matrix is None:
            M = matrix(self._module.base_ring(), len(self._points),
                       self._module.degree(), self._points)
            M.set_immutable()
            self._matrix = M
        return self._matrix

    def module(self):
        r"""
        Return the ambient module of ``self``.

        OUTPUT:

        - a :class:`module <FreeModule_generic>`.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.module()
            3-d lattice N
        """
        return self._module

    @classmethod    # @staticmethod does not work in Cython so far
    def output_format(cls, format=None):
        r"""
        Return or set the output format for **ALL** point collections.

        INPUT:

        - ``format`` -- (optional) if given, must be one of the strings
            * "default" -- output one point per line with vertical alignment of
              coordinates in text mode, same as "tuple" for LaTeX;
            * "tuple" -- output ``tuple(self)`` with lattice information;
            * "matrix" -- output :meth:`matrix` with lattice information;
            * "column matrix" -- output :meth:`column_matrix` with lattice
              information;
            * "separated column matrix" -- same as "column matrix" for text
              mode, for LaTeX separate columns by lines (not shown by jsMath).

        OUTPUT:

        - a string with the current format (only if ``format`` was omitted).

        This function affects both regular and LaTeX output.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c
            N(0, 0, 1),
            N(1, 0, 1),
            N(0, 1, 1),
            N(1, 1, 1)
            in 3-d lattice N
            sage: c.output_format()
            'default'
            sage: c.output_format("tuple")
            sage: c
            (N(0, 0, 1), N(1, 0, 1), N(0, 1, 1), N(1, 1, 1))
            in 3-d lattice N
            sage: c.output_format("matrix")
            sage: c
            [0 0 1]
            [1 0 1]
            [0 1 1]
            [1 1 1]
            in 3-d lattice N
            sage: c.output_format("column matrix")
            sage: c
            [0 1 0 1]
            [0 0 1 1]
            [1 1 1 1]
            in 3-d lattice N
            sage: c.output_format("separated column matrix")
            sage: c
            [0 1 0 1]
            [0 0 1 1]
            [1 1 1 1]
            in 3-d lattice N

        Note that the last two outpus are identical, separators are only
        inserted in the LaTeX mode::

            sage: latex(c)
            \left(\begin{array}{r|r|r|r}
            0 & 1 & 0 & 1 \\
            0 & 0 & 1 & 1 \\
            1 & 1 & 1 & 1
            \end{array}\right)_{N}

        Since this is a static method, you can call it for the class directly::

            sage: from sage.geometry.point_collection import PointCollection
            sage: PointCollection.output_format("default")
            sage: c
            N(0, 0, 1),
            N(1, 0, 1),
            N(0, 1, 1),
            N(1, 1, 1)
            in 3-d lattice N
        """
        global _output_format
        if format is None:
            return _output_format
        assert format in ["default", "tuple", "matrix", "column matrix",
                          "separated column matrix"]
        _output_format = format

    def set(self):
        r"""
        Return points of ``self`` as a :class:`frozenset`.

        OUTPUT:

        - a :class:`frozenset`.

        EXAMPLES::

            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)]).rays()
            sage: c.set()
            frozenset({N(0, 0, 1), N(0, 1, 1), N(1, 0, 1), N(1, 1, 1)})
        """
        if self._set is None:
            self._set = frozenset(self._points)
        return self._set
