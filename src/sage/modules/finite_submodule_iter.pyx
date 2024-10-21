# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.rings.finite_rings
r"""
Iterators over finite submodules of a `\ZZ`-module

We iterate over the elements of a finite `\ZZ`-module. The action
of `\ZZ` must be the natural one.

This class is intended to provide optimizations for the
:meth:`sage.free_module.FreeModule_generic:__iter__` method.

AUTHORS:

- Thomas Feulner (2012-08-31): initial version
- Punarbasu Purkayastha (2012-11-09): replaced the loop with recursion
- Thomas Feulner (2012-11-09): added functionality to enumerate cosets, FiniteFieldsubspace_projPoint_iterator

EXAMPLES::

    sage: from sage.modules.finite_submodule_iter import FiniteZZsubmodule_iterator
    sage: F.<x,y,z> = FreeAlgebra(GF(3),3)
    sage: iter = FiniteZZsubmodule_iterator([x,y], [3,3])
    sage: list(iter)
    [0, x, 2*x, y, x + y, 2*x + y, 2*y, x + 2*y, 2*x + 2*y]

There is a specialization for subspaces over finite fields::

    sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_iterator
    sage: A = random_matrix(GF(4, 'a'), 5, 100)
    sage: iter = FiniteFieldsubspace_iterator(A)
    sage: len(list(iter))
    1024

The module also allows the iteration over cosets::

    sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_iterator
    sage: A = random_matrix(GF(4, 'a'), 5, 100)
    sage: v = random_vector(GF(4, 'a'), 100)
    sage: iter = FiniteFieldsubspace_iterator(A, v)
    sage: len(list(iter))
    1024

TESTS::

    sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_iterator
    sage: A = random_matrix(GF(4, 'a'), 5, 100)
    sage: iter = FiniteFieldsubspace_iterator(A)
    sage: TestSuite(iter).run(skip='_test_pickling')

In a second step, we will replace all calls to ``__iter__`` for finite submodules. This
will result in improved running times::

    sage: A = random_matrix(GF(2), 15, 100)
    sage: X = A.row_space()
    sage: x = [0 for _ in X] #long time #takes 7.12 seconds
    sage: y = [0 for _ in FiniteFieldsubspace_iterator(A)] # takes 0.05 seconds
    sage: sorted(x) == sorted(y) #long time
    True
"""

#*****************************************************************************
#       Copyright (C) 2012 Thomas Feulner <thomas.feulner@uni-bayreuth.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef class FiniteZZsubmodule_iterator:
    r"""
    Let `G` be an abelian group and suppose that `(g_0, \ldots, g_n)`
    is a list of elements of `G`, whose additive orders are equal to `m_i`
    and `\sum_{i=0}^n x_i g_i = 0` for `x_i \in \ZZ_{m_i}` for
    `i \in \{0, \dots, n\}` implies `x_i=0` for all `i`.

    This class implements an iterator over the `\ZZ`-submodule
    `M = \{\sum_{i=0}^n x_i g_i\}`. If the independence condition from
    above is not fulfilled, we can still use this iterator to run over the
    elements. In this case the elements will occur multiple times.

    Getting from one element of the submodule to another is performed by
    one single addition in `G`.

    INPUT:

        - ``basis`` -- the elements `(g_0, \ldots, g_n)`
        - ``order`` -- (optional) the additive_orders `m_i` of `g_i`
        - ``coset_rep`` -- (optional) an element of `g`,
          if one aims to compute a coset of the `\ZZ`-submodule `M`
        - ``immutable`` -- boolean (default: ``False``); set it to
          ``True`` to return immutable elements. Setting this to
          ``True`` makes sense if the elements are vectors. See
          :class:`FiniteFieldsubspace_iterator` for examples.

    EXAMPLES::

        sage: from sage.modules.finite_submodule_iter import FiniteZZsubmodule_iterator
        sage: F.<x,y,z> = FreeAlgebra(GF(3),3)
        sage: iter = FiniteZZsubmodule_iterator([x,y], [3,3])
        sage: list(iter)
        [0, x, 2*x, y, x + y, 2*x + y, 2*y, x + 2*y, 2*x + 2*y]
        sage: iter = FiniteZZsubmodule_iterator([x,y], [3,3], z)
        sage: list(iter)
        [z, x + z, 2*x + z, y + z, x + y + z, 2*x + y + z, 2*y + z, x + 2*y + z, 2*x + 2*y + z]
    """

    def __init__(self, basis, order=None, coset_rep=None, immutable=False):
        """
        See :class:`FiniteZZsubmodule_iterator`.

        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteZZsubmodule_iterator
            sage: F.<x,y,z> = FreeAlgebra(GF(3),3)
            sage: iter = FiniteZZsubmodule_iterator([x,y], [3,3])
            sage: list(iter)
            [0, x, 2*x, y, x + y, 2*x + y, 2*y, x + 2*y, 2*x + 2*y]
        """
        if order is None:
            try:
                order = [b.additive_order() for b in basis]
            except (AttributeError, NotImplementedError):
                raise ValueError("Unable to determine the additive order "
                                 "of a basis element. Use the optional "
                                 "parameter `order`.")

        self._basis = basis[0]
        self._basis_all = basis
        self._basis_length = len(basis)
        self._count = 0
        self._immutable = immutable

        if coset_rep is None:
            self._coset_rep = self._basis.parent().zero()
        else:
            self._coset_rep = self._basis.parent()(coset_rep)
        if self._basis_length == 1:
            self._cw = self._coset_rep
        else:
            self._cw = self._basis.parent().zero()
            self._other_ZZ = FiniteZZsubmodule_iterator(basis[1:], order[1:], coset_rep)

        self._order = order[0]
        self._other = self._basis.parent().zero()  # dummy initialization
        self._plus = [self._cw]  # storing this provides about 20% speedup
        for _ in range(self._order - 1):
            self._cw += self._basis
            self._plus.append(self._cw)

    def __next__(self):
        """
        Return the next submodule element. This will just add/subtract
        another element of the ``basis``.

        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteZZsubmodule_iterator
            sage: F.<x,y,z> = FreeAlgebra(GF(3),3)
            sage: iter = FiniteZZsubmodule_iterator([x,y], [3,3])
            sage: next(iter) #indirect doctest
            0
            sage: next(iter) #indirect doctest
            x
        """
        self._iteration()
        v = self._cw
        if self._immutable:
            v.set_immutable()
        return v

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteZZsubmodule_iterator
            sage: F.<x,y,z> = FreeAlgebra(GF(3),3)
            sage: FiniteZZsubmodule_iterator([x,y], [3,3])
            Iterator over ZZ-submodule generated by [x, y]
        """
        return "Iterator over ZZ-submodule generated by {0}".format(self._basis_all)

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteZZsubmodule_iterator
            sage: F.<x,y,z> = FreeAlgebra(GF(3),3)
            sage: list(FiniteZZsubmodule_iterator([x,y], [3,3])) #indirect doctest
            [0, x, 2*x, y, x + y, 2*x + y, 2*y, x + 2*y, 2*x + 2*y]
        """
        return self

    cdef ModuleElement _iteration(FiniteZZsubmodule_iterator self):
        """
        This is the core implementation of the iteration.

        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteZZsubmodule_iterator
            sage: F.<x,y,z> = FreeAlgebra(GF(3),3)
            sage: iter = FiniteZZsubmodule_iterator([x,y], [3,3])
            sage: next(iter) #indirect doctest
            0
            sage: next(iter), next(iter), next(iter) #indirect doctest
            (x, 2*x, y)
        """
        if self._basis_length == 1:
            if self._count < self._order:
                self._cw = self._plus[self._count]
                self._count += 1
            else:
                raise StopIteration
        else:
            if self._count == 0 or self._count == self._order:
                self._other = next(self._other_ZZ)
                self._cw = < ModuleElement > self._other
                self._count = 1
            else:
                self._cw = self._other + self._plus[self._count]
                self._count += 1


cdef class FiniteFieldsubspace_iterator(FiniteZZsubmodule_iterator):
    """
    This class implements an iterator over the subspace of a vector
    space over a finite field. The subspace is generated by ``basis``.

    INPUT:

        - ``basis`` -- list of vectors or a matrix with elements over
          a finite field. If a matrix is provided then it is not checked
          whether the matrix is full ranked. Similarly, if a list of
          vectors is provided, then the linear independence of the vectors
          is not checked.

        - ``coset_rep`` -- (optional) a vector in the same ambient space,
          if one aims to compute a coset of the vector space given by ``basis``

        - ``immutable`` -- boolean (default: ``False``); set it to
          ``True`` to return immutable vectors

    EXAMPLES::

        sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_iterator
        sage: A = random_matrix(GF(2), 10, 100)
        sage: iter = FiniteFieldsubspace_iterator(A)
        sage: len(list(iter))
        1024
        sage: X = random_matrix(GF(4, 'a'), 7, 100).row_space()
        sage: s = list(X)  # long time (5s on sage.math, 2013)
        sage: t = list(FiniteFieldsubspace_iterator(X.basis()))  # takes 0.31s
        sage: sorted(t) == sorted(s)  # long time
        True

    TESTS:

    We test whether we get immutable vectors if immutable=True::

        sage: iter = FiniteFieldsubspace_iterator(A, immutable=True)
        sage: c = next(iter)
        sage: c.is_immutable()
        True
    """

    def __init__(self, basis, coset_rep=None, immutable=False):
        """
        See :class:`FiniteFieldsubspace_iterator`.

        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_iterator
            sage: A = random_matrix(GF(2), 10, 100)
            sage: iter = FiniteFieldsubspace_iterator(A)
            sage: X = list(iter)
            sage: len(X)
            1024
            sage: v = random_vector(GF(2), 100)
            sage: iter = FiniteFieldsubspace_iterator(A, v)
            sage: Y = list(iter)
            sage: len(Y)
            1024
            sage: all(Y[i]-X[i] == v for i in range(len(X)))
            True
        """
        cdef Py_ssize_t i, p
        cdef list pows, order

        F = basis[0].base_ring()
        P = F.prime_subfield()
        p = P.order()
        alpha = F.primitive_element()
        degree = F.degree()

        pows = [alpha ** i for i in range(degree)]
        basis = [_p * x for x in basis for _p in pows]  # a ZZ_p-basis for the vectorspace
        order = [p] * (len(basis))

        FiniteZZsubmodule_iterator.__init__(self, basis, order, coset_rep,
                                            immutable=immutable)

cdef class FiniteFieldsubspace_projPoint_iterator:
    """
    This class implements an iterator over the projective points of a vector
    space over a finite field. The vector space is generated by ``basis`` and
    need not to be equal to the full ambient space.

    A projective point (= one dimensional subspace) `P` will be represented by a
    generator `p`. To ensure that all `p` will be normalized you can set the
    optional argument ``normalize`` to ``True``.

    INPUT:

        - ``basis`` -- list of vectors or a matrix with elements over
          a finite field. If a matrix is provided then it is not checked
          whether the matrix is full ranked. Similarly, if a list of
          vectors is provided, then the linear independence of the vectors
          is not checked.

        - ``normalize`` -- boolean (default: ``False``); whether the returned
          vectors should be normalized, i.e. the first nonzero coordinate is
          equal to 1

        - ``immutable`` -- boolean (default: ``False``); set it to
          ``True`` to return immutable vectors

    EXAMPLES::

        sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_iterator, FiniteFieldsubspace_projPoint_iterator
        sage: A = random_matrix(GF(4, 'a'), 5, 100)
        sage: a = len(list(FiniteFieldsubspace_iterator(A)))
        sage: b = len(list(FiniteFieldsubspace_projPoint_iterator(A)))
        sage: b == (a-1)/3
        True

    Prove that the option ``normalize == True`` will only return normalized vectors.

        sage: all(x.monic() == x for x in FiniteFieldsubspace_projPoint_iterator(A, True))
        True

    TESTS::

        sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_projPoint_iterator
        sage: A = MatrixSpace(GF(7), 10, 10).one()
        sage: len(list(FiniteFieldsubspace_projPoint_iterator(A[:0]))) #indirect doctest
        0
        sage: len(list(FiniteFieldsubspace_projPoint_iterator(A[:1]))) #indirect doctest
        1
        sage: len(list(FiniteFieldsubspace_projPoint_iterator(A[:2]))) #indirect doctest
        8
        sage: iter = FiniteFieldsubspace_projPoint_iterator(A[:2], immutable=True)
        sage: next(iter).is_immutable()
        True
    """

    def __init__(self, basis, normalize=False, immutable=False):
        """
        See :class:`FiniteFieldsubspace_projPoint_iterator`.

        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_projPoint_iterator
            sage: A = random_matrix(GF(4, 'a'), 4, 100)
            sage: iter = FiniteFieldsubspace_projPoint_iterator(A)
            sage: len(list(iter))
            85
        """
        from sage.matrix.constructor import matrix
        self._basis = list(basis)
        self._basis_length = len(self._basis)
        self._immutable = immutable
        if immutable:
            for b in self._basis:
                b.set_immutable()
        if normalize:
            B = matrix(self._basis)
            B.echelonize()
            self._basis = B.rows()
            self._basis.reverse()

        if self._basis_length == 0:
            self._one_dimensional_case = 2
        else:
            self._one_dimensional_case = 1

    def __next__(self):
        """
        Return the next projective point. This will just add/subtract
        another element of the ``basis`` except for the cases when the rank will
        increase.

        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_projPoint_iterator
            sage: A = MatrixSpace(GF(3), 10,10).one()
            sage: iter = FiniteFieldsubspace_projPoint_iterator(A)
            sage: next(iter) #indirect doctest
            (1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: next(iter) #indirect doctest
            (0, 1, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        if self._one_dimensional_case > 0:
            if self._one_dimensional_case == 1:
                self._one_dimensional_case = 2
                # this returns immutable vectors if immutable is True
                return self._basis[0]
            else:
                if self._basis_length > 1:
                    self._it = FiniteFieldsubspace_iterator(self._basis[:1],
                                                            self._basis[1],
                                                            immutable=self._immutable)
                    self._normalized_pos = 1
                    self._one_dimensional_case = 0
                else:
                    raise StopIteration
        try:
            return next(self._it)
        except StopIteration:
            self._normalized_pos += 1
            if self._normalized_pos == self._basis_length:
                raise StopIteration
            else:
                self._it = FiniteFieldsubspace_iterator(self._basis[:self._normalized_pos],
                                                        self._basis[self._normalized_pos],
                                                        immutable=self._immutable)
                return next(self._it)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_projPoint_iterator
            sage: A = MatrixSpace(GF(3), 2, 2).one()
            sage: FiniteFieldsubspace_projPoint_iterator(A)
            Iterator over the projective points of a subspace generated by [(1, 0), (0, 1)]
        """

        return "Iterator over the projective points of a subspace generated by {0}".format(self._basis)

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.modules.finite_submodule_iter import FiniteFieldsubspace_projPoint_iterator
            sage: A = MatrixSpace(GF(3), 10,10).one()
            sage: len(list(FiniteFieldsubspace_projPoint_iterator(A))) #indirect doctest
            29524
            sage: A = MatrixSpace(GF(3), 1,1).one()
            sage: len(list(FiniteFieldsubspace_projPoint_iterator(A))) #indirect doctest
            1
        """
        return self
