r"""
Super Lie Algebras

These are `K`-graded Lie algebras satisfying graded antisymmetry and
the graded Jacobi identity, where there is a natural map `K \to \ZZ / 2 \ZZ`.

Note that these are more commonly referred to as Lie superalgebras in the
literature.

AUTHORS:

- Aditya Dwarkesh, Martin Frankland (08-08-2023): Initial version
- Travis Scrimshaw (25-02-2024): Refactoring with preparation for future
  implementations
"""

# ****************************************************************************
#       Copyright (C) 2023 Aditya Dwarkesh <ad19ms047@iiserkol.ac.in>
#                          Martin Frankland <Martin.Frankland@uregina.ca>
#                     2024 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.indexed_generators import standardize_names_index_set
from sage.algebras.lie_algebras.structure_coefficients import WithStructureCoefficients_abstract, _standardize_s_coeff
from sage.algebras.lie_algebras.lie_algebra_element import SuperStructureCoefficientsElement
from sage.sets.family import Family
from sage.categories.lie_algebras import LieAlgebras


class SuperLieAlgebra(Parent, UniqueRepresentation):
    r"""
    A super Lie algebra.

    Currently this only supports finite dimensional Lie superalgebras
    specified by their structure coefficients with integer degrees for
    each basis element.

    INPUT:

    - ``R`` -- the base ring
    - ``s_coeff`` -- the structure coefficients
    - ``names`` -- (optional) the names of the basis elements
    - ``index_set`` -- (optional) the indices of the basis elements
    - ``degrees`` -- (default: all 0) list of degrees of the basis elements

    EXAMPLES::

        sage: d = {('x','y'): {'z':1}}
        sage: L.<x,y,z> = SuperLieAlgebra(QQ, d, degrees=(1,1,2))  # corresponding to Q^2 + Q
    """
    @staticmethod
    def __classcall_private__(cls, R=None, arg0=None, arg1=None, degrees=None, names=None,
                              index_set=None, category=None, **kwds):
        """
        Determine the correct subclass to construct based on the input.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, (1,1,2))
            sage: type(L)
            <class 'sage.algebras.lie_algebras.superliealgebra.SuperLieAlgebraStructureCoefficients_with_category'>
        """
        # Currently only defined by structure coefficients is implemented.
        # Therefore, we assume the input is valid for that.
        if degrees is None:
            degrees = arg1
            arg1 = None
        return SuperLieAlgebraStructureCoefficients(R, arg0, degrees=degrees, names=names, index_set=index_set, **kwds)

    def __getitem__(self, x):
        """
        If ``x`` is a pair `(a, b)`, return the Lie bracket `[a, b]
        (including if `a` or `b` are Lie (sub)algebras, in which case the
        corresponding ideal is constructed).
        Otherwise try to return the `x`-th element of ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
            sage: L[x, y]
            z
        """
        if isinstance(x, tuple) and len(x) == 2:
            # Check if we need to construct an ideal
            Cat = LieAlgebras.Super
            if x[0] in Cat:
                if x[1] in Cat:
                    #return x[0].product_space(x[1])
                    raise NotImplementedError("sub Lie superalgebras not yet implemented")
                #return x[0].ideal(x[1])
                raise NotImplementedError("ideals not yet implemented")
            elif x[1] in Cat:
                #return x[1].ideal(x[0])
                raise NotImplementedError("ideals not yet implemented")
            # Otherwise it is the bracket of two elements
            return self(x[0])._bracket_(self(x[1]))
        return super().__getitem__(x)

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
            sage: elt = L([x, y]); elt
            z
            sage: elt.parent() is L
            True

        TESTS:

        Check that `0` gives the zero element::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
            sage: L(0)
            0
            sage: L(0).parent() is L
            True
        """
        if isinstance(x, list) and len(x) == 2:
            return self(x[0])._bracket_(self(x[1]))

        if x == 0:
            return self.zero()

        try:
            if x in self.module():
                return self.from_vector(x)
        except AttributeError:
            pass

        if x in self.base_ring():
            # We have already handled the case when x == 0
            raise ValueError("can only convert the scalar 0 into a Lie algebra element")

        return self.element_class(self, x)

    def _Hom_(self, Y, category):
        """
        Return the homset from ``self`` to ``Y`` in the category ``category``.

        INPUT:

        - ``Y`` -- a Lie superalgebra
        - ``category`` -- a subcategory of :class:`LieAlgebras` or ``None``

        The sole purpose of this method is to construct the homset
        as a :class:`~sage.algebras.lie_algebras.morphism.LieAlgebraHomset`.

        This method is not meant to be called directly. Please use
        :func:`sage.categories.homset.Hom` instead.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, (1,1,2))
            sage: Hom(L, L)
            Set of Super Lie algebra morphisms from
             Super Lie algebra generated by (x, y, z) in degrees (1, 1, 2) over Rational Field
             to Super Lie algebra generated by (x, y, z) in degrees (1, 1, 2) over Rational Field
        """
        cat = LieAlgebras(self.base_ring()).Super()
        if category is not None and not category.is_subcategory(cat):
            raise TypeError(f"{category} is not a subcategory of Lie superalgebras")
        if Y not in cat:
            raise TypeError(f"{Y} is not a Lie superalgebra")
        from sage.algebras.lie_algebras.morphism import LieAlgebraHomset
        return LieAlgebraHomset(self, Y, category=category)


class SuperLieAlgebraStructureCoefficients(WithStructureCoefficients_abstract, SuperLieAlgebra):
    """
    A Lie superalgebra defined in terms of structure coefficients.
    """
    @staticmethod
    def __classcall_private__(cls, R, s_coeff, degrees=None, names=None, index_set=None, **kwargs):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: L1.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, [1,1,2])
            sage: L2 = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, [1,1,2], names=['x','y','z'])
            sage: L1 is L2
            True
        """
        if names is None:
            if degrees is None:
                raise ValueError("you must specify names or degrees")
            else:
                n = len(degrees)
            names = tuple('x{}'.format(i) for i in range(n))
        elif isinstance(names, str):
            names = tuple(names.split(','))
            n = len(names)
        else:
            n = len(names)
            names = tuple(names)

        if degrees is None:
            degrees = tuple([0] * n)

        degrees = tuple(degrees)
        names, index_set = standardize_names_index_set(names, index_set)
        s_coeff = _standardize_s_coeff(s_coeff, index_set, degrees)

        return super().__classcall__(cls, R, s_coeff, degrees, names, index_set, **kwargs)

    def __init__(self, R, s_coeff, degrees, names, index_set, category=None, **kwds):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, (1,1,2))
            sage: TestSuite(L).run()
        """
        self._degrees = degrees
        WithStructureCoefficients_abstract.__init__(self, R, s_coeff, names, index_set, **kwds)
        category = LieAlgebras(R).WithBasis().FiniteDimensional().Super().or_subcategory(category)
        SuperLieAlgebra.__init__(self, base=R, names=names, category=category)
        self.__ngens = len(self._indices)

    def _test_degree(self, **options):
        """
        Test to ensure that the bracket satisfies the grading property.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, (1,1,2))
            sage: L._test_degree()

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z':1}}, (1,1,4))
            sage: L._test_degree()
            Traceback (most recent call last):
            ...
            AssertionError: 2 != 4
        """
        tester = self._tester(**options)
        dim = len(self._indices)
        # TODO: Verify this at initialization
        for i in range(dim):
            for j in range(i, dim):
                k = (i, j)
                if (i, j) in self._s_coeff:
                    elt = self.from_vector(self._s_coeff[k])
                    tester.assertTrue(elt.is_homogeneous())
                    tester.assertEqual(self._degrees[i] + self._degrees[j],
                                       elt.degree())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
            sage: L
            Super Lie algebra generated by (x, y, z) in degrees (1, 1, 2) over Rational Field
        """
        return "Super Lie algebra generated by {} in degrees {} over {}".format(self.gens(), self._degrees, self.base_ring())

    def degree_on_basis(self, x):
        """
        Return the degree of the basis element indexed by ``x``.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
            sage: L.indices()
            {'x', 'y', 'z'}
            sage: [L.degree_on_basis(i) for i in L.indices()]
            [1, 1, 2]
        """
        return self._degrees[self._index_to_pos[x]]

    @cached_method
    def super_lie_algebra_generators(self):
        """
        Return the generators of ``self`` as a super Lie algebra.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
            sage: L.super_lie_algebra_generators()
            Finite family {'x': x, 'y': y, 'z': z}
        """
        return Family(self._indices, self.monomial, name="monomial map")

    @cached_method
    def gens(self):
        r"""
        Return a tuple whose entries are the generators for this
        object, in some order.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
            sage: L.gens()
            (x, y, z)
        """
        G = self.super_lie_algebra_generators()
        try:
            return tuple(G[i] for i in self.variable_names())
        except (KeyError, IndexError):
            return tuple(G[i] for i in self.indices())
        except ValueError:
            return tuple(G)

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = SuperLieAlgebra(QQ, {('x','y'): {'z': 1}}, (1,1,2))
            sage: L.gen(0)
            x
        """
        return self.gens()[i]

    class Element(SuperStructureCoefficientsElement, WithStructureCoefficients_abstract.Element):
        pass
