r"""
Toric lattice elements

This module was designed as a part of the framework for toric varieties
(:mod:`~sage.schemes.toric.variety`,
:mod:`~sage.schemes.toric.fano_variety`).

AUTHORS:

- Andrey Novoseltsev (2010-05-27): initial version.

TESTS:

Let's create some lattices first::

    sage: N = ToricLattice(3)
    sage: M = N.dual()

Now we are ready to create elements of toric lattices::

    sage: n = N([1,2,3])
    sage: n
    N(1, 2, 3)
    sage: m = M(1,2,3)
    sage: m
    M(1, 2, 3)

Dual lattices can act on each other::

    sage: n * m
    14
    sage: m * n
    14

You can also add elements of the same lattice or scale them::

    sage: 2 * n
    N(2, 4, 6)
    sage: n * 2
    N(2, 4, 6)
    sage: n + n
    N(2, 4, 6)

However, you cannot "mix wrong lattices" in your expressions::

    sage: n + m
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for +:
    '3-d lattice N' and '3-d lattice M'
    sage: n * n
    Traceback (most recent call last):
    ...
    TypeError: elements of the same toric lattice cannot be multiplied!
    sage: n == m
    False

Note that ``n`` and ``m`` are not equal to each other even though they are
both "just (1,2,3)." Moreover, you cannot easily convert elements between
toric lattices::

    sage: M(n)
    Traceback (most recent call last):
    ...
    TypeError: N(1, 2, 3) cannot be converted to 3-d lattice M!

If you really need to consider elements of one lattice as elements of another,
you can either use intermediate conversion to "just a vector"::

    sage: ZZ3 = ZZ^3
    sage: n_in_M = M(ZZ3(n))
    sage: n_in_M
    M(1, 2, 3)
    sage: n == n_in_M
    False
    sage: n_in_M == m
    True

Or you can create a homomorphism from one lattice to any other::

    sage: h = N.hom(identity_matrix(3), M)
    sage: h(n)
    M(1, 2, 3)
"""
# The "tutorial" above is a truncated version of one in toric_lattice.py.


#*****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.mpz cimport *

from sage.geometry.toric_plotter import ToricPlotter
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.structure.coerce_exceptions import CoercionException
from sage.structure.element cimport Vector
from sage.rings.integer cimport Integer
from sage.structure.richcmp cimport richcmp_not_equal


def is_ToricLatticeElement(x):
    r"""
    Check if ``x`` is an element of a toric lattice.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is an element of a toric lattice, ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.toric_lattice_element import (
        ....:   is_ToricLatticeElement)
        sage: is_ToricLatticeElement(1)
        False
        sage: e = ToricLattice(3).an_element()
        sage: e
        N(1, 0, 0)
        sage: is_ToricLatticeElement(e)
        True
    """
    return isinstance(x, ToricLatticeElement)


# Why do we need a special class:
# - customize output to include lattice name
# - prohibit operations mixing "wrong" lattices
cdef class ToricLatticeElement(Vector_integer_dense):
    r"""
    Create an element of a toric lattice.

    .. WARNING::

        You probably should not construct such elements explicitly.

    INPUT:

    - same as for
      :class:`~sage.modules.vector_integer_dense.Vector_integer_dense`.

    OUTPUT:

    - element of a toric lattice.

    TESTS::

        sage: N = ToricLattice(3)
        sage: from sage.geometry.toric_lattice_element import (
        ....:           ToricLatticeElement)
        sage: e = ToricLatticeElement(N, [1,2,3])
        sage: e
        N(1, 2, 3)
        sage: TestSuite(e).run()
    """
    def __richcmp__(self, right, op):
        r"""
        Compare ``self`` and ``right`` according to the operator ``op``.

        INPUT:

        - ``right`` -- another ToricLatticeElement

        OUTPUT:

        boolean

        First compare the ambient toric lattice, then compare the vectors.

        TESTS::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: n = N(1,2,3)
            sage: m = M(1,2,3)
            sage: n == m
            False
            sage: n2 = N(1,2,3)
            sage: n == n2
            True
            sage: n is n2
            False
        """
        if not is_ToricLatticeElement(right):
            return NotImplemented

        PL_ambient = self.parent().ambient_module()
        PR_ambient = right.parent().ambient_module()
        if PL_ambient != PR_ambient:
            return richcmp_not_equal(PL_ambient, PR_ambient, op)
        # Now use the real comparison of vectors
        return self._richcmp_(right, op)

    # For some reason, vectors work just fine without redefining this function
    # from the base class, but if it is not here, we get "unhashable type"...
    def __hash__(self):
        r"""
        Return the hash of ``self``.

        OUTPUT:

        - integer.

        TESTS::

            sage: N = ToricLattice(3)
            sage: n = N(1,2,3)
            sage: hash(n)
            Traceback (most recent call last):
            ...
            TypeError: mutable vectors are unhashable
            sage: n.set_immutable()
            sage: hash(n) == hash(n)
            True
        """
        return Vector_integer_dense.__hash__(self)

    cpdef _act_on_(self, other, bint self_on_left):
        """
        Act on ``other``.

        INPUT:

        - ``other`` - :class:`ToricLatticeElement`.

        OUTPUT:

        - integer, if ``other`` is an element of the dual lattice of ``self``;

        - ``CoercionException`` is raised if ``other`` is an element of
          an incompatible toric lattice;

        - standard output for ``self`` acting as an integral vector on
          ``other`` if the latter one is not an element of a toric lattice.

        TESTS::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: n = N(1,2,3)
            sage: m = M(1,2,3)
            sage: n * m # indirect doctest
            14

        Now we test behaviour with other types::

            sage: v = vector([1, 2, 3])
            sage: v * n == n * v
            True
            sage: v = vector([1, 1/2, 3/4])
            sage: v * n == n * v
            True
            sage: A = matrix(3, range(9))
            sage: A * n
            (8, 26, 44)
            sage: n * A
            (24, 30, 36)
            sage: B = A / 3
            sage: B * n
            (8/3, 26/3, 44/3)
            sage: n * B
            (8, 10, 12)
        """
        Ns = self.parent()
        # We try to deal only with the case of two lattice elements...
        if is_ToricLatticeElement(other):
            if other.parent().ambient_module() is Ns.ambient_module().dual():
                # Our own _dot_product_ is disabled
                return Vector_integer_dense._dot_product_(self, other)
            raise CoercionException("only elements of dual toric lattices "
                                    "can act on each other!")
        # ... however we also need to treat the case when other is an integral
        # vector, since otherwise it will be coerced to the parent of self and
        # then the dot product will be called for elements of the same lattice
        if isinstance(other, Vector_integer_dense):
            return Vector_integer_dense._dot_product_(self, other)
        # We also allow action on elements of lattice quotients
        try:
            lift = other.lift()
            if is_ToricLatticeElement(lift):
                if other.parent().W().is_submodule(Ns.dual().W()):
                    return Vector_integer_dense._dot_product_(self, lift)
                raise CoercionException("only elements of dual toric lattices "
                                        "can act on each other!")
        except AttributeError:  # No lift
            pass
        # Now let the standard framework work...
        return Vector_integer_dense._act_on_(self, other, self_on_left)

    # We need to override this function to prohibit default behaviour.
    # It seems to be called when right is in the same lattice as self, which
    # is wrong from our point of view.
    cpdef _dot_product_(self, Vector right):
        """
        Raise a ``TypeError`` exception.

        Dot product is not defined on toric lattices (there are actions of
        dual lattices on each other instead).

        INPUT:

        - ``right`` - vector.

        OUTPUT:

        - ``TypeError`` exception is raised.

        TESTS::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: n = N(1,2,3)
            sage: m = M(1,2,3)
            sage: n * n # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: elements of the same
            toric lattice cannot be multiplied!
        """
        raise TypeError("elements of the same toric lattice cannot be "
                        "multiplied!")

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: Ld = ToricLattice(3, "L").dual()
            sage: e = Ld(1,2,3)
            sage: e._latex_()
            '\\left(1,\\,2,\\,3\\right)_{L^*}'
        """
        return "%s_{%s}" % (super()._latex_(),
                            self.parent().ambient_module()._latex_name)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: Ld = ToricLattice(3, "L").dual()
            sage: e = Ld(1,2,3)
            sage: e._repr_()
            'L*(1, 2, 3)'
        """
        return (self.parent().ambient_module()._name
                + super()._repr_())

    def __reduce__(self):
        """
        Override the base ``__reduce__`` to correctly pickle/unpickle elements.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: loads(dumps(N(1,2,3)))
            N(1, 2, 3)
        """
        return (unpickle_v1, (self._parent, self.list(), self._degree,
                              not self._is_immutable))

    def plot(self, **options):
        r"""
        Plot ``self``.

        INPUT:

        - any options for toric plots (see :func:`toric_plotter.options
          <sage.geometry.toric_plotter.options>`), none are mandatory.

        OUTPUT:

        - a plot.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: n = N(1,2,3)
            sage: n.plot()                                                              # needs sage.plot
            Graphics3d Object
        """
        tp = ToricPlotter(options, self.parent().degree())
        tp.adjust_options()
        return tp.plot_points([self])


def unpickle_v1(parent, entries, degree, is_mutable):
    """
    Unpickle a :class:`ToricLatticeElement`

    INPUT:

    - ``parent`` -- The parent toric lattice.

    - ``entries`` -- a list. The coordinates of the lattice point.

    - ``degree`` -- integer. the dimension of the toric lattice.

    - ``is_mutable`` -- boolean. Whether the lattice element is mutable.

    OUTPUT:

    The :class:`ToricLatticeElement` determined by the input data.

    EXAMPLES::

        sage: N = ToricLattice(3, "lattice")
        sage: loads(dumps(N(1,2,3)))   # indirect test
        lattice(1, 2, 3)
        sage: from sage.geometry.toric_lattice_element import unpickle_v1
        sage: unpickle_v1(N,[1,2,3],3,True)
        lattice(1, 2, 3)
    """
    cdef ToricLatticeElement v
    v = ToricLatticeElement.__new__(ToricLatticeElement)
    v._init(degree, parent)
    cdef Integer z
    cdef Py_ssize_t i
    for i in range(degree):
        z = Integer(entries[i])
        mpz_set(v._entries[i], z.value)
    v._is_immutable = not is_mutable
    return v
