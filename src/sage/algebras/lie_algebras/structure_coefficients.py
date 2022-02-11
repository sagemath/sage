"""
Lie Algebras Given By Structure Coefficients

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
#from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.indexed_generators import (IndexedGenerators,
                                               standardize_names_index_set)

from sage.categories.lie_algebras import LieAlgebras

from sage.algebras.lie_algebras.lie_algebra_element import StructureCoefficientsElement
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra
#from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
#from sage.algebras.lie_algebras.ideal import LieAlgebraIdeal
#from sage.algebras.lie_algebras.quotient import QuotientLieAlgebra
from sage.modules.free_module import FreeModule
from sage.sets.family import Family

class LieAlgebraWithStructureCoefficients(FinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    A Lie algebra with a set of specified structure coefficients.

    The structure coefficients are specified as a dictionary `d` whose
    keys are pairs of basis indices, and whose values are
    dictionaries which in turn are indexed by basis indices. The value
    of `d` at a pair `(u, v)` of basis indices is the dictionary whose
    `w`-th entry (for `w` a basis index) is the coefficient of `b_w`
    in the Lie bracket `[b_u, b_v]` (where `b_x` means the basis
    element with index `x`).

    INPUT:

    - ``R`` -- a ring, to be used as the base ring

    - ``s_coeff`` -- a dictionary, indexed by pairs of basis indices
      (see below), and whose values are dictionaries which are
      indexed by (single) basis indices and whose values are elements
      of `R`

    - ``names`` -- list or tuple of strings

    - ``index_set`` -- (default: ``names``) list or tuple of hashable
      and comparable elements

    OUTPUT:

    A Lie algebra over ``R`` which (as an `R`-module) is free with
    a basis indexed by the elements of ``index_set``. The `i`-th
    basis element is displayed using the name ``names[i]``.
    If we let `b_i` denote this `i`-th basis element, then the Lie
    bracket is given by the requirement that the `b_k`-coefficient
    of `[b_i, b_j]` is ``s_coeff[(i, j)][k]`` if
    ``s_coeff[(i, j)]`` exists, otherwise ``-s_coeff[(j, i)][k]``
    if ``s_coeff[(j, i)]`` exists, otherwise `0`.

    EXAMPLES:

    We create the Lie algebra of `\QQ^3` under the Lie bracket defined
    by `\times` (cross-product)::

        sage: L = LieAlgebra(QQ, 'x,y,z', {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}})
        sage: (x,y,z) = L.gens()
        sage: L.bracket(x, y)
        z
        sage: L.bracket(y, x)
        -z

    TESTS:

    We can input structure coefficients that fail the Jacobi
    identity, but the test suite will call us out on it::

        sage: Fake = LieAlgebra(QQ, 'x,y,z', {('x','y'):{'z':3}, ('y','z'):{'z':1}, ('z','x'):{'y':1}})
        sage: TestSuite(Fake).run()
        Failure in _test_jacobi_identity:
        ...

    Old tests !!!!!placeholder for now!!!!!::

        sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
        sage: L.basis()
        Finite family {'x': x, 'y': y}
    """
    @staticmethod
    def __classcall_private__(cls, R, s_coeff, names=None, index_set=None, **kwds):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: L1 = LieAlgebra(QQ, 'x,y', {('x','y'): {'x':1}})
            sage: L2 = LieAlgebra(QQ, 'x,y', {('y','x'): {'x':-1}})
            sage: L1 is L2
            True

        Check that we convert names to the indexing set::

            sage: L = LieAlgebra(QQ, 'x,y,z', {('x','y'): {'z':1}, ('y','z'): {'x':1}, ('z','x'): {'y':1}}, index_set=list(range(3)))
            sage: (x,y,z) = L.gens()
            sage: L[x,y]
            L[2]
        """
        names, index_set = standardize_names_index_set(names, index_set)

        # Make sure the structure coefficients are given by the index set
        if names is not None and names != tuple(index_set):
            d = {x: index_set[i] for i,x in enumerate(names)}
            get_pairs = lambda X: X.items() if isinstance(X, dict) else X
            try:
                s_coeff = {(d[k[0]], d[k[1]]): [(d[x], y) for x,y in get_pairs(s_coeff[k])]
                           for k in s_coeff}
            except (KeyError, ValueError):
                # At this point we assume they are given by the index set
                pass

        s_coeff = LieAlgebraWithStructureCoefficients._standardize_s_coeff(s_coeff, index_set)
        if s_coeff.cardinality() == 0:
            from sage.algebras.lie_algebras.abelian import AbelianLieAlgebra
            return AbelianLieAlgebra(R, names, index_set, **kwds)

        if (names is None and len(index_set) <= 1) or len(names) <= 1:
            from sage.algebras.lie_algebras.abelian import AbelianLieAlgebra
            return AbelianLieAlgebra(R, names, index_set, **kwds)

        return super(LieAlgebraWithStructureCoefficients, cls).__classcall__(
            cls, R, s_coeff, names, index_set, **kwds)

    @staticmethod
    def _standardize_s_coeff(s_coeff, index_set):
        """
        Helper function to standardize ``s_coeff`` into the appropriate form
        (dictionary indexed by pairs, whose values are dictionaries).
        Strips items with coefficients of 0 and duplicate entries.
        This does not check the Jacobi relation (nor antisymmetry if the
        cardinality is infinite).

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
            sage: d = {('y','x'): {'x':-1}}
            sage: LieAlgebraWithStructureCoefficients._standardize_s_coeff(d, ('x', 'y'))
            Finite family {('x', 'y'): (('x', 1),)}
        """
        # Try to handle infinite basis (once/if supported)
        #if isinstance(s_coeff, AbstractFamily) and s_coeff.cardinality() == infinity:
        #    return s_coeff

        index_to_pos = {k: i for i,k in enumerate(index_set)}

        sc = {}
        # Make sure the first gen is smaller than the second in each key
        for k in s_coeff.keys():
            v = s_coeff[k]
            if isinstance(v, dict):
                v = v.items()

            if index_to_pos[k[0]] > index_to_pos[k[1]]:
                key = (k[1], k[0])
                vals = tuple((g, -val) for g, val in v if val != 0)
            else:
                if not index_to_pos[k[0]] < index_to_pos[k[1]]:
                    if k[0] == k[1]:
                        if not all(val == 0 for g, val in v):
                            raise ValueError("elements {} are equal but their bracket is not set to 0".format(k))
                        continue
                key = tuple(k)
                vals = tuple((g, val) for g, val in v if val != 0)

            if key in sc.keys() and sorted(sc[key]) != sorted(vals):
                raise ValueError("two distinct values given for one and the same bracket")

            if vals:
                sc[key] = vals
        return Family(sc)

    def __init__(self, R, s_coeff, names, index_set, category=None, prefix=None,
                 bracket=None, latex_bracket=None, string_quotes=None, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'): {'x':1}})
            sage: TestSuite(L).run()
        """
        default = (names != tuple(index_set))
        if prefix is None:
            if default:
                prefix = 'L'
            else:
                prefix = ''
        if bracket is None:
            bracket = default
        if latex_bracket is None:
            latex_bracket = default
        if string_quotes is None:
            string_quotes = default

        #self._pos_to_index = dict(enumerate(index_set))
        self._index_to_pos = {k: i for i,k in enumerate(index_set)}
        if "sorting_key" not in kwds:
            kwds["sorting_key"] = self._index_to_pos.__getitem__

        cat = LieAlgebras(R).WithBasis().FiniteDimensional().or_subcategory(category)
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, index_set, cat)
        IndexedGenerators.__init__(self, self._indices, prefix=prefix,
                                   bracket=bracket, latex_bracket=latex_bracket,
                                   string_quotes=string_quotes, **kwds)

        self._M = FreeModule(R, len(index_set))

        # Transform the values in the structure coefficients to elements
        def to_vector(tuples):
            vec = [R.zero()]*len(index_set)
            for k,c in tuples:
                vec[self._index_to_pos[k]] = c
            vec = self._M(vec)
            vec.set_immutable()
            return vec
        self._s_coeff = {(self._index_to_pos[k[0]], self._index_to_pos[k[1]]):
                         to_vector(s_coeff[k])
                         for k in s_coeff.keys()}

    # For compatibility with CombinatorialFreeModuleElement
    _repr_term = IndexedGenerators._repr_generator
    _latex_term = IndexedGenerators._latex_generator

    def structure_coefficients(self, include_zeros=False):
        """
        Return the dictionary of structure coefficients of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z', {('x','y'): {'x':1}})
            sage: L.structure_coefficients()
            Finite family {('x', 'y'): x}
            sage: S = L.structure_coefficients(True); S
            Finite family {('x', 'y'): x, ('x', 'z'): 0, ('y', 'z'): 0}
            sage: S['x','z'].parent() is L
            True

        TESTS:

        Check that :trac:`23373` is fixed::

            sage: L = lie_algebras.sl(QQ, 2)
            sage: sorted(L.structure_coefficients(True), key=str)
            [-2*E[-alpha[1]], -2*E[alpha[1]], h1]
        """
        if not include_zeros:
            pos_to_index = dict(enumerate(self._indices))
            return Family({(pos_to_index[k[0]], pos_to_index[k[1]]):
                           self.element_class(self, self._s_coeff[k])
                           for k in self._s_coeff})
        ret = {}
        zero = self._M.zero()
        for i,x in enumerate(self._indices):
            for j, y in enumerate(self._indices[i+1:]):
                if (i, j+i+1) in self._s_coeff:
                    elt = self._s_coeff[i, j+i+1]
                elif (j+i+1, i) in self._s_coeff:
                    elt = -self._s_coeff[j+i+1, i]
                else:
                    elt = zero
                ret[x,y] = self.element_class(self, elt) # +i+1 for offset
        return Family(ret)

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.dimension()
            2
        """
        return self.basis().cardinality()

    def module(self, sparse=True):
        """
        Return ``self`` as a free module.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}})
            sage: L.module()
            Sparse vector space of dimension 3 over Rational Field
        """
        return FreeModule(self.base_ring(), self.dimension(), sparse=sparse)

    @cached_method
    def zero(self):
        """
        Return the element `0` in ``self``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: L.zero()
            0
        """
        return self.element_class(self, self._M.zero())

    def monomial(self, k):
        """
        Return the monomial indexed by ``k``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: L.monomial('x')
            x
        """
        return self.element_class(self, self._M.basis()[self._index_to_pos[k]])

    def term(self, k, c=None):
        """
        Return the term indexed by ``i`` with coefficient ``c``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: L.term('x', 4)
            4*x
        """
        if c is None:
            c = self.base_ring().one()
        else:
            c = self.base_ring()(c)
        return self.element_class(self, c * self._M.basis()[self._index_to_pos[k]])

    def from_vector(self, v):
        """
        Return an element of ``self`` from the vector ``v``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
            sage: L.from_vector([1, 2, -2])
            x + 2*y - 2*z
        """
        return self.element_class(self, self._M(v))

    def some_elements(self):
        """
        Return some elements of ``self``.

        EXAMPLES::

            sage: L = lie_algebras.three_dimensional(QQ, 4, 1, -1, 2)
            sage: L.some_elements()
            [X, Y, Z, X + Y + Z]
        """
        return list(self.basis()) + [self.sum(self.basis())]

    def change_ring(self, R):
        r"""
        Return a Lie algebra with identical structure coefficients over ``R``.

        INPUT:

        - ``R`` -- a ring

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(ZZ, {('x','y'): {'z':1}})
            sage: L.structure_coefficients()
            Finite family {('x', 'y'): z}
            sage: LQQ = L.change_ring(QQ)
            sage: LQQ.structure_coefficients()
            Finite family {('x', 'y'): z}
            sage: LSR = LQQ.change_ring(SR)
            sage: LSR.structure_coefficients()
            Finite family {('x', 'y'): z}
        """
        return LieAlgebraWithStructureCoefficients(
            R, self.structure_coefficients(),
            names=self.variable_names(), index_set=self.indices())

    class Element(StructureCoefficientsElement):
        def _sorted_items_for_printing(self):
            """
            Return a list of pairs ``(k, c)`` used in printing.

            .. WARNING::

                The internal representation order is fixed, whereas this
                depends on ``"sorting_key"`` print option as it is used
                only for printing.

            EXAMPLES::

                sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'): {'z':1}})
                sage: elt = x + y/2 - z; elt
                x + 1/2*y - z
                sage: elt._sorted_items_for_printing()
                [('x', 1), ('y', 1/2), ('z', -1)]
                sage: key = {'x': 2, 'y': 1, 'z': 0}
                sage: L.print_options(sorting_key=key.__getitem__)
                sage: elt._sorted_items_for_printing()
                [('z', -1), ('y', 1/2), ('x', 1)]
                sage: elt
                -z + 1/2*y + x
            """
            print_options = self.parent().print_options()
            pos_to_index = dict(enumerate(self.parent()._indices))
            v = [(pos_to_index[k], c) for k, c in self.value.items()]
            try:
                v.sort(key=lambda monomial_coeff:
                            print_options['sorting_key'](monomial_coeff[0]),
                       reverse=print_options['sorting_reverse'])
            except Exception: # Sorting the output is a plus, but if we can't, no big deal
                pass
            return v

