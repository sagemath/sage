.. -*- coding: utf-8 -*-
.. _tutorial-implementing-algebraic-structures:

===========================================
Tutorial: Implementing Algebraic Structures
===========================================

.. linkall

.. MODULEAUTHOR:: Nicolas M. Thiéry <nthiery at users.sf.net>,
   Jason Bandlow <jbandlow@gmail.com> et al.

This tutorial will cover four concepts:

* endowing free modules and vector spaces with additional algebraic structure
* defining morphisms
* defining coercions and conversions
* implementing algebraic structures with several realizations

At the end of this tutorial, the reader should be able to reimplement
by himself the example of algebra with several realizations::

    sage: Sets().WithRealizations().example()
    The subset algebra of {1, 2, 3} over Rational Field

Namely, we consider an algebra `A(S)` whose basis is indexed by the
subsets `s` of a given set `S`. `A(S)` is endowed with three natural
basis: ``F``, ``In``, ``Out``; in the first basis, the product is
given by the union of the indexing sets. The ``In`` basis and ``Out``
basis are defined respectively by:

    .. MATH::

        In_s  = \sum_{t\subset s} F_t \qquad
        F_s   = \sum_{t\supset s} Out_t

Each such basis gives a realization of `A`, where the elements are
represented by their expansion in this basis. In the running exercises
we will progressively implement this algebra and its three
realizations, with coercions and mixed arithmetic between them.

This tutorial heavily depends on :ref:`sage.modules.tutorial_free_modules`.
You may also want to read the less specialized thematic tutorial
:ref:`How to implement new algebraic structures <coercion_and_categories>`.

Subclassing free modules and including category information
===========================================================

As a warm-up, we implement the group algebra of the additive group
`\ZZ/5\ZZ`. Of course this is solely for pedagogical purposes; group
algebras are already implemented (see ``ZMod(5).algebra(ZZ)``). Recall
that a fully functional `\ZZ`-module over this group can be created
with the simple command::

    sage: A = CombinatorialFreeModule(ZZ, Zmod(5), prefix='a')

We reproduce the same, but by deriving a subclass of
:class:`CombinatorialFreeModule`::

    sage: class MyCyclicGroupModule(CombinatorialFreeModule):
    ....:     """An absolutely minimal implementation of a module whose basis is a cyclic group"""
    ....:     def __init__(self, R, n, *args, **kwargs):
    ....:         CombinatorialFreeModule.__init__(self, R, Zmod(n), *args, **kwargs)

    sage: A = MyCyclicGroupModule(QQ, 6, prefix='a') # or 4 or 5 or 11     ...
    sage: a = A.basis()
    sage: A.an_element()
    2*a[0] + 2*a[1] + 3*a[2]

We now want to endow `A` with its natural product structure, to get
the desired group algebra. To define a multiplication, we should be in
a category where multiplication makes sense, which is not yet the
case::

    sage: A.category()
    Category of finite dimensional vector spaces with basis over Rational Field

We can look at the available `Categories <../reference/categories/sage/categories/category.html>`_ 
from the documentation in the reference manual or we can use introspection to
look through the list of categories to pick one we want::

    sage: sage.categories.<tab>                   # not tested

Once we have chosen an appropriate category (here
:class:`AlgebrasWithBasis`), one can look at one example::

    sage: E = AlgebrasWithBasis(QQ).example(); E
    An example of an algebra with basis: the free algebra on the generators ('a', 'b', 'c') over Rational Field
    sage: e = E.an_element(); e
    B[word: ] + 2*B[word: a] + 3*B[word: b] + B[word: bab]

and browse through its code::

    sage: E??                                     # not tested

This code is meant as a template for implementing a
new algebra. In particular, this template suggests that we need to implement the
methods ``product_on_basis``, ``one_basis``, ``_repr_`` and
``algebra_generators``. Another way to get this list of methods is to
ask the category (TODO: find a slicker idiom for this)::

    sage: from sage.misc.abstract_method import abstract_methods_of_class
    sage: abstract_methods_of_class(AlgebrasWithBasis(QQ).element_class)
    {'optional': ['_add_', '_mul_'],
     'required': ['__bool__', 'monomial_coefficients']}
    sage: abstract_methods_of_class(AlgebrasWithBasis(QQ).parent_class)
    {'optional': ['one_basis', 'product_on_basis'], 'required': ['__contains__']}

.. WARNING::

    The result above is not yet necessarily complete; many required
    methods in the categories are not yet marked as
    :func:`abstract_methods`. We also recommend browsing the
    documentation of this category: :class:`AlgebrasWithBasis`.

Adding these methods, here is the minimal implementation of the group algebra::

    sage: class MyCyclicGroupAlgebra(CombinatorialFreeModule):
    ....:
    ....:     def __init__(self, R, n, **keywords):
    ....:         self._group = Zmod(n)
    ....:         CombinatorialFreeModule.__init__(self, R, self._group,
    ....:             category=AlgebrasWithBasis(R), **keywords)
    ....:
    ....:     def product_on_basis(self, left, right):
    ....:         return self.monomial( left + right )
    ....:
    ....:     def one_basis(self):
    ....:         return self._group.zero()
    ....:
    ....:     def algebra_generators(self):
    ....:         return Family( [self.monomial( self._group(1) ) ] )
    ....:
    ....:     def _repr_(self):
    ....:         return "Jason's group algebra of %s over %s"%(self._group, self.base_ring())

Some notes about this implementation:

* Alternatively, we could have defined ``product`` instead of
  ``product_on_basis``::

   ....:     # def product(self, left, right):
   ....:     #     return ## something ##

* For the sake of readability in this tutorial, we have stripped out
  all the documentation strings. Of course all of those should be
  present as in ``E``.

* The purpose of ``**keywords`` is to pass down options like
  ``prefix`` to :class:`CombinatorialFreeModules`.


Let us do some calculations::

    sage: A = MyCyclicGroupAlgebra(QQ, 2, prefix='a') # or 4 or 5 or 11     ...
    sage: a = A.basis();
    sage: f = A.an_element();
    sage: A, f
    (Jason's group algebra of Ring of integers modulo 2 over Rational Field, 2*a[0] + 2*a[1])
    sage: f * f
    8*a[0] + 8*a[1]
    sage: f.<tab>                                 # not tested
    sage: f.is_idempotent()
    False
    sage: A.one()
    a[0]
    sage: x = A.algebra_generators().first() # Typically x,y,    ... = A.algebra_generators()
    sage: [x^i for i in range(4)]
    [a[0], a[1], a[0], a[1]]
    sage: g = 2*a[1]; (f + g)*f == f*f + g*f
    True

This seems to work fine, but we would like to put more stress on our
implementation to shake potential bugs out of it. To this end, we will
use :class:`TestSuite`, a tool that performs many routine tests
on our algebra for us.

Since we defined the class interactively, instead of in a Python
module, those tests will complain about "pickling". We can silence this
error by making sage think that the class is defined in a module. We could also
just ignore those failing tests for now or call :class:`TestSuite` with the
argument ``skip='_test_pickling')``::

    sage: import __main__
    sage: __main__.MyCyclicGroupAlgebra = MyCyclicGroupAlgebra

Ok, let's run the tests::

    sage: TestSuite(A).run(verbose=True)
    running ._test_additive_associativity() . . . pass
    running ._test_an_element() . . . pass
    running ._test_associativity() . . . pass
    running ._test_cardinality() . . . pass
    running ._test_category() . . . pass
    running ._test_characteristic() . . . pass
    running ._test_construction() . . . pass
    running ._test_distributivity() . . . pass
    running ._test_elements() . . .
      Running the test suite of self.an_element()
      running ._test_category() . . . pass
      running ._test_eq() . . . pass
      running ._test_monomial_coefficients() . . . pass
      running ._test_new() . . . pass
      running ._test_nonzero_equal() . . . pass
      running ._test_not_implemented_methods() . . . pass
      running ._test_pickling() . . . pass
      pass
    running ._test_elements_eq_reflexive() . . . pass
    running ._test_elements_eq_symmetric() . . . pass
    running ._test_elements_eq_transitive() . . . pass
    running ._test_elements_neq() . . . pass
    running ._test_eq() . . . pass
    running ._test_new() . . . pass
    running ._test_not_implemented_methods() . . . pass
    running ._test_one() . . . pass
    running ._test_pickling() . . . pass
    running ._test_prod() . . . pass
    running ._test_some_elements() . . . pass
    running ._test_zero() . . . pass

For more information on categories, see :ref:`sage.categories.primer`::

    sage: sage.categories.primer?                 # not tested

Review
------

We wanted to implement an algebra, so we:

#.  Created the underlying vector space using :class:`CombinatorialFreeModule`
#.  Looked at ``sage.categories.<tab>`` to find an appropriate category
#.  Loaded an example of that category, and used :func:`sage.misc.abstract_method.abstract_methods_of_class`, to see what methods we needed to write
#.  Added the category information and other necessary methods to our class
#.  Ran :class:`TestSuite` to catch potential discrepancies

Exercises
---------

#. Make a tiny modification to ``product_on_basis`` in
   "MyCyclicGroupAlgebra" to implement the *dual* of the group algebra
   of the cyclic group instead of its group algebra (so the product is now given by
   `b_fb_g=\delta_{f,g}b_f`).

   Run the :class:`TestSuite` tests (you may ignore the "pickling"
   errors). What do you notice?

   Fix the implementation of ``one`` and check that the :class:`TestSuite` tests now pass.

   Add the Hopf algebra structure. Hint: look at the example::

       sage: C = HopfAlgebrasWithBasis(QQ).example()


#. Given a set `S`, say::

        sage: S = Set([1,2,3,4,5])

   and a base ring, say::

        sage: R = QQ

   implement an `R`-algebra::

        sage: F = SubsetAlgebraOnFundamentalBasis(S, R)   # todo: not implemented

   with a basis ``(b_s)_{s\subset S}`` indexed by the subsets of ``S``::

        sage: Subsets(S)
        Subsets of {1, 2, 3, 4, 5}

   and where the product is defined by `b_s b_t = b_{s\cup t}`.


Morphisms
=========

To better understand relationships between algebraic spaces, one wants
to consider morphisms between them::

    sage: A.module_morphism?                      # not tested
    sage: A = MyCyclicGroupAlgebra(QQ, 2, prefix='a')
    sage: B = MyCyclicGroupAlgebra(QQ, 6, prefix='b')
    sage: A, B
    (Jason's group algebra of Ring of integers modulo 2 over Rational Field, Jason's group algebra of Ring of integers modulo 6 over Rational Field)

::

    sage: def func_on_basis(g):
    ....:     r"""
    ....:     This function is the 'brain' of a (linear) morphism
    ....:     from A --> B.
    ....:     The input is the index of basis element of the domain (A).
    ....:     The output is an element of the codomain (B).
    ....:     """
    ....:     if g==1: return B.monomial(Zmod(6)(3))# g==1 in the range A
    ....:     else:    return B.one()

We can now define a morphism that extends this function to `A` by
linearity::

    sage: phi = A.module_morphism(func_on_basis, codomain=B)
    sage: f = A.an_element()
    sage: f
    2*a[0] + 2*a[1]
    sage: phi(f)
    2*b[0] + 2*b[3]


Exercise
--------

Define a new free module ``In`` with basis indexed by the subsets of
`S`, and a morphism ``phi`` from ``In`` to ``F`` defined by

    .. MATH:: \phi(In_s) = \sum_{t\subset s} F_t


Diagonal and Triangular Morphisms
=================================

We now illustrate how to specify that a given morphism is diagonal or triangular
with respect to some order on the basis, which means that the morphism is
invertible and Sage is able to compute the inverse morphism automatically.
Currently this feature requires the domain and codomain to have the same index
set (in progress ...).

::

    sage: X = CombinatorialFreeModule(QQ, Partitions(), prefix='x'); x = X.basis();
    sage: Y = CombinatorialFreeModule(QQ, Partitions(), prefix='y'); y = Y.basis();

A diagonal module morphism takes as argument a function whose input is
the index of a basis element of the domain, and whose output is the
coefficient of the corresponding basis element of the codomain::

    sage: def diag_func(p):
    ....:     if len(p)==0: return 1
    ....:     else: return p[0]
    ....:
    ....:
    sage: diag_func(Partition([3,2,1]))
    3
    sage: X_to_Y = X.module_morphism(diagonal=diag_func, codomain=Y)
    sage: f = X.an_element();
    sage: f
    2*x[[]] + 2*x[[1]] + 3*x[[2]]
    sage: X_to_Y(f)
    2*y[[]] + 2*y[[1]] + 6*y[[2]]

Python fun fact: ``~`` is the inversion operator (but be careful with
int's!)::

    sage: ~2
    1/2
    sage: ~(int(2)) # in python this is the bitwise complement: ~x = -x-1
    -3

Diagonal module morphisms are invertible::

    sage: Y_to_X = ~X_to_Y
    sage: f = y[Partition([3])] - 2*y[Partition([2,1])]
    sage: f
    -2*y[[2, 1]] + y[[3]]
    sage: Y_to_X(f)
    -x[[2, 1]] + 1/3*x[[3]]
    sage: X_to_Y(Y_to_X(f))
    -2*y[[2, 1]] + y[[3]]

For triangular morphisms, just like ordinary morphisms, we need a
function that accepts as input the index of a basis element of the
domain and returns an element of the codomain.  We think of this
function as representing the columns of the matrix of the linear
transformation::

    sage: def triang_on_basis(p):
    ....:     return Y.sum_of_monomials(mu for mu in Partitions(sum(p)) if mu >= p)
    ....:
    sage: triang_on_basis([3,2])
    y[[3, 2]] + y[[4, 1]] + y[[5]]
    sage: X_to_Y = X.module_morphism(triang_on_basis, triangular='lower', unitriangular=True, codomain=Y)
    sage: f = x[Partition([1,1,1])] + 2*x[Partition([3,2])];
    sage: f
    x[[1, 1, 1]] + 2*x[[3, 2]]

::

    sage: X_to_Y(f)
    y[[1, 1, 1]] + y[[2, 1]] + y[[3]] + 2*y[[3, 2]] + 2*y[[4, 1]] + 2*y[[5]]

Triangular module_morphisms are also invertible, even if ``X`` and
``Y`` are both infinite-dimensional::

    sage: Y_to_X = ~X_to_Y
    sage: f
    x[[1, 1, 1]] + 2*x[[3, 2]]
    sage: Y_to_X(X_to_Y(f))
    x[[1, 1, 1]] + 2*x[[3, 2]]

For details, see
:meth:`ModulesWithBasis.ParentMethods.module_morphism` (and also
:class:`sage.categories.modules_with_basis.TriangularModuleMorphism`)::

    sage: A.module_morphism?                      # not tested

Exercise
--------

Redefine the morphism ``phi`` from the previous exercise as a morphism that is
triangular with respect to inclusion of subsets and define the inverse morphism.
You may want to use the following comparison key as
``key`` argument to ``modules_morphism``::

    sage: def subset_key(s):
    ....:     """
    ....:     A comparison key on sets that gives a linear extension
    ....:     of the inclusion order.
    ....:
    ....:     INPUT:
    ....:
    ....:      - ``s`` -- set
    ....:
    ....:     EXAMPLES::
    ....:
    ....:         sage: sorted(Subsets([1,2,3]), key=subset_key)
    ....:         [{}, {1}, {2}, {3}, {1, 2}, {1, 3}, {2, 3}, {1, 2, 3}]
    ....:     """
    ....:     return (len(s), list(s))


Coercions
=========

Once we have defined a morphism from `X \to Y`, we can register it as
a coercion.  This will allow Sage to apply the morphism automatically
whenever we combine elements of `X` and `Y` together. See
http://sagemath.com/doc/reference/coercion.html for more
information. As a training step, let us first define a morphism `X` to
`Y`, and register it as a coercion::

    sage: def triang_on_basis(p):
    ....:     return Y.sum_of_monomials(mu for mu in Partitions(sum(p)) if mu >= p)

    sage: triang_on_basis([3,2])
    y[[3, 2]] + y[[4, 1]] + y[[5]]
    sage: X_to_Y = X.module_morphism(triang_on_basis, triangular='lower', unitriangular=True, codomain=Y)
    sage: X_to_Y.<tab>                            # not tested
    sage: X_to_Y.register_as_coercion()

Now we can not only convert elements from `X` to `Y`, but we can also do
mixed arithmetic with these elements::

    sage: Y(x[Partition([3,2])])
    y[[3, 2]] + y[[4, 1]] + y[[5]]
    sage: Y([2,2,1]) + x[Partition([2,2,1])]
    2*y[[2, 2, 1]] + y[[3, 1, 1]] + y[[3, 2]] + y[[4, 1]] + y[[5]]


Exercise
--------

Use the inverse of ``phi`` to implement the inverse coercion from
``F`` to ``In``. Reimplement ``In`` as an algebra, with a product
method making it use ``phi`` and its inverse.


A digression: new bases and quotients of symmetric functions
============================================================

As an application, we show how to combine what we have learned to
implement a new basis and a quotient of the algebra of symmetric
functions::

    sage: SF = SymmetricFunctions(QQ);  # A graded Hopf algebra
    sage: h  = SF.homogeneous()         # A particular basis, indexed by partitions (with some additional magic)

So, `h` is a graded algebra whose basis is indexed by partitions. In more
detail, ``h([i])`` is the sum of all monomials of degree `i`::

    sage: h([2]).expand(4)
    x0^2 + x0*x1 + x1^2 + x0*x2 + x1*x2 + x2^2 + x0*x3 + x1*x3 + x2*x3 + x3^2

and ``h(mu) = prod( h(p) for p in mu )``::

    sage: h([3,2,2,1]) == h([3]) * h([2]) * h([2]) * h([1])
    True

Here we define a new basis `(X_\lambda)_\lambda` by triangularity
with respect to `h`; namely, we set `X_\lambda = \sum_{\mu\geq \lambda, |\mu|=|\nu|} h_\mu`::

    sage: class MySFBasis(CombinatorialFreeModule):
    ....:     r"""
    ....:     Note: We would typically use SymmetricFunctionAlgebra_generic
    ....:     for this. This is as an example only.
    ....:     """
    ....:
    ....:     def __init__(self, R, *args, **kwargs):
    ....:         """ TODO: Informative doc-string and examples """
    ....:         CombinatorialFreeModule.__init__(self, R, Partitions(), category=AlgebrasWithBasis(R), *args, **kwargs)
    ....:         self._h = SymmetricFunctions(R).homogeneous()
    ....:         self._to_h = self.module_morphism( self._to_h_on_basis, triangular='lower', unitriangular=True, codomain=self._h)
    ....:         self._from_h = ~(self._to_h)
    ....:         self._to_h.register_as_coercion()
    ....:         self._from_h.register_as_coercion()
    ....:
    ....:     def _to_h_on_basis(self, la):
    ....:         return self._h.sum_of_monomials(mu for mu in Partitions(sum(la)) if mu >= la)
    ....:
    ....:     def product(self, left, right):
    ....:         return self( self._h(left) * self._h(right) )
    ....:
    ....:     def _repr_(self):
    ....:         return "Jason's basis for symmetric functions over %s"%self.base_ring()
    ....:
    ....:     @cached_method
    ....:     def one_basis(self):
    ....:         r""" Returns the index of the basis element that is equal to '1'."""
    ....:         return Partition([])
    sage: X = MySFBasis(QQ, prefix='x'); x = X.basis(); h = SymmetricFunctions(QQ).homogeneous()
    sage: f = X(h([2,1,1]) - 2*h([2,2]))  # Note the capital X
    sage: f
    x[[2, 1, 1]] - 3*x[[2, 2]] + 2*x[[3, 1]]
    sage: h(f)
    h[2, 1, 1] - 2*h[2, 2]
    sage: f*f*f
    x[[2, 2, 2, 1, 1, 1, 1, 1, 1]] - 7*x[[2, 2, 2, 2, 1, 1, 1, 1]] + 18*x[[2, 2, 2, 2, 2, 1, 1]]
    - 20*x[[2, 2, 2, 2, 2, 2]] + 8*x[[3, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
    sage: h(f*f)
    h[2, 2, 1, 1, 1, 1] - 4*h[2, 2, 2, 1, 1] + 4*h[2, 2, 2, 2]

We now implement a quotient of the algebra of symmetric functions
obtained by killing any monomial symmetric function `m_\lambda` such
that the first part of `\lambda` is greater than `k`. See
:meth:`Sets.SubcategoryMethods.Subquotients` for more details about
implementing quotients::

    sage: class MySFQuotient(CombinatorialFreeModule):
    ....:     r"""
    ....:     The quotient of the ring of symmetric functions by the ideal generated
    ....:     by those monomial symmetric functions whose part is larger than some fixed
    ....:     number ``k``.
    ....:     """
    ....:     def __init__(self, R, k, prefix=None, *args, **kwargs):
    ....:         CombinatorialFreeModule.__init__(self, R,
    ....:             Partitions(NonNegativeIntegers(), max_part=k),
    ....:             prefix = 'mm',
    ....:             category = Algebras(R).Graded().WithBasis().Quotients(), *args, **kwargs)
    ....:
    ....:         self._k = k
    ....:         self._m = SymmetricFunctions(R).monomial()
    ....:
    ....:         self.lift = self.module_morphism(self._m.monomial)
    ....:         self.retract = self._m.module_morphism(self._retract_on_basis, codomain=self)
    ....:
    ....:         self.lift.register_as_coercion()
    ....:         self.retract.register_as_coercion()
    ....:
    ....:     def ambient(self):
    ....:         return self._m
    ....:
    ....:     def _retract_on_basis(self, mu):
    ....:         r"""
    ....:         Takes the index of a basis element of a monomial
    ....:         symmetric function, and returns the projection of that
    ....:         element to the quotient.
    ....:         """
    ....:         if len(mu) > 0 and mu[0] > self._k:
    ....:             return self.zero()
    ....:         return self.monomial(mu)
    ....:
    sage: MM = MySFQuotient(QQ, 3)
    sage: mm = MM.basis()
    sage: m = SymmetricFunctions(QQ).monomial()
    sage: P = Partition
    sage: g = m[P([3,2,1])] + 2*m[P([3,3])] + m[P([4,2])]; g
    m[3, 2, 1] + 2*m[3, 3] + m[4, 2]
    sage: f = MM(g); f
    mm[[3, 2, 1]] + 2*mm[[3, 3]]
    sage: m(f)
    m[3, 2, 1] + 2*m[3, 3]

    sage: (m(f))^2
    8*m[3, 3, 2, 2, 1, 1] + 12*m[3, 3, 2, 2, 2] + 24*m[3, 3, 3, 2, 1] + 48*m[3, 3, 3, 3]
    + 4*m[4, 3, 2, 2, 1] + 4*m[4, 3, 3, 1, 1] + 14*m[4, 3, 3, 2] + 4*m[4, 4, 2, 2]
    + 4*m[4, 4, 3, 1] + 6*m[4, 4, 4] + 4*m[5, 3, 2, 1, 1] + 4*m[5, 3, 2, 2]
    + 12*m[5, 3, 3, 1] + 2*m[5, 4, 2, 1] + 6*m[5, 4, 3] + 4*m[5, 5, 1, 1] + 2*m[5, 5, 2]
    + 4*m[6, 2, 2, 1, 1] + 6*m[6, 2, 2, 2] + 6*m[6, 3, 2, 1] + 10*m[6, 3, 3] + 2*m[6, 4, 1, 1] + 5*m[6, 4, 2] + 4*m[6, 5, 1] + 4*m[6, 6]

    sage: f^2
    8*mm[[3, 3, 2, 2, 1, 1]] + 12*mm[[3, 3, 2, 2, 2]] + 24*mm[[3, 3, 3, 2, 1]] + 48*mm[[3, 3, 3, 3]]

    sage: (m(f))^2 - m(f^2)
    4*m[4, 3, 2, 2, 1] + 4*m[4, 3, 3, 1, 1] + 14*m[4, 3, 3, 2] + 4*m[4, 4, 2, 2] + 4*m[4, 4, 3, 1] + 6*m[4, 4, 4] + 4*m[5, 3, 2, 1, 1] + 4*m[5, 3, 2, 2] + 12*m[5, 3, 3, 1] + 2*m[5, 4, 2, 1] + 6*m[5, 4, 3] + 4*m[5, 5, 1, 1] + 2*m[5, 5, 2] + 4*m[6, 2, 2, 1, 1] + 6*m[6, 2, 2, 2] + 6*m[6, 3, 2, 1] + 10*m[6, 3, 3] + 2*m[6, 4, 1, 1] + 5*m[6, 4, 2] + 4*m[6, 5, 1] + 4*m[6, 6]

    sage: MM( (m(f))^2 - m(f^2) )
    0

Implementing algebraic structures with several realizations
===========================================================

We now return to the subset algebra and use it as an example to show how to
implement several different bases for an algebra with automatic coercions
between the different bases. We have already implemented three bases for this
algebra:  the ``F``, ``In``, and ``Out`` bases, as well as coercions between
them. In real calculations it is convenient to tie these parents together by
implementing an object ``A`` that models the abstract algebra itself. Then, the
parents ``F``, ``In`` and ``Out`` will be *realizations* of ``A``, while ``A``
will be a *parent with realizations*. See :func:`Sets().WithRealizations
<sage.categories.with_realizations.WithRealizations>` for more information
about the expected user interface and the rationale.

Here is a brief template highlighting the overall structure:

.. CODE-BLOCK:: python

    class MyAlgebra(Parent, UniqueRepresentation):
        def __init__(self, R, ...):
            category = Algebras(R).Commutative()
            Parent.__init__(self, category=category.WithRealizations())
            # attribute initialization, construction of the morphisms
            # between the bases, ...

        class Bases(Category_realization_of_parent):
            def super_categories(self):
                A = self.base()
                category = Algebras(A.base_ring()).Commutative()
                return [A.Realizations(), category.Realizations().WithBasis()]

            class ParentMethods:
                r"""Code that is common to all bases of the algebra"""

            class ElementMethods:
                r"""Code that is common to elements of all bases of the algebra"""

        class FirstBasis(CombinatorialFreeModule, BindableClass):
            def __init__(self, A):
                CombinatorialFreeModule.__init__(self, ..., category=A.Bases())

            # implementation of the multiplication, the unit, ...

        class SecondBasis(CombinatorialFreeModule, BindableClass):
            def __init__(self, A):
                CombinatorialFreeModule.__init__(self, ..., category=A.Bases())

            # implementation of the multiplication, the unit, ...


The class ``MyAlgebra`` implements a commutative algebra ``A`` with several
realizations, which we specify in the constructor of ``MyAlgebra``. The two
bases classes ``MyAlgebra.FirstBasis`` and ``MyAlgebra.SecondBasis`` implement
different realizations of ``A`` that correspond to distinguished bases on which
elements are expanded. They are initialized in the category ``MyAlgebra.Bases``
of all bases of ``A``, whose role is to factor out their common features. In
particular, this construction says that they are:

- realizations of ``A``
- realizations of a commutative algebra, with a distinguished basis

.. NOTE::

    There is a bit of redundancy here: given that ``A`` knows it is a
    commutative algebra with realizations the infrastructure could, in
    principle, determine that its realizations are commutative algebras. If this
    was done then it would be possible to implement ``Bases.super_categories`` by
    returning::

            [A.Realizations().WithBasis()]

    However, this has not been implemented yet.

.. NOTE::

    Inheriting from :class:`BindableCass` just provides syntactic
    sugar: it makes ``MyAlgebras().FirstBasis()`` a shorthand for
    ``MyAlgebras.FirstBasis(MyAlgebras().FirstBasis())`` (binding
    behavior). The class ``Bases`` inherits this binding behavior from
    :class:`Category_realization_of_parent` , which is why we can
    write ``MyAlgebras().Bases`` instead of
    ``MyAlgebras.Bases(MyAlgebras())``

.. NOTE::

    More often than not, the constructors for all of the bases will be very
    similar, if not identical; so we would want to factor it out. Annoyingly,
    the natural approach of putting the constructor in ``Bases.ParentMethods``
    does not work because this is an abstract class whereas the constructor
    handles the concrete implementation of the data structure. Similarly, it
    would be better if it was only necessary to  specify the classes the bases
    inherit from once, but this can't code go into ``Bases`` for the same
    reason.

    The current recommended solution is to have an additional class ``Basis``
    that factors out the common concrete features of the different bases:

    .. CODE-BLOCK:: python

        ...

        class Basis(CombinatorialFreeModule, BindableClass):
            def __init__(self, A):
                CombinatorialFreeModule.__init__(self, ..., category=A.Bases())

        class FirstBasis(Basis):
            ...

        class SecondBasis(Basis):
            ...

    This solution works but it is not optimal because to share features between
    the two bases code needs to go into two locations, ``Basis`` and ``Bases``,
    depending on whether they are concrete or abstract, respectively.

We now urge the reader to browse the full code of the following
example, which is meant as a complete template for constructing new
parents with realizations::

    sage: A = Sets().WithRealizations().example(); A
    The subset algebra of {1, 2, 3} over Rational Field

    sage: A??                                     # not implemented


Review
======

Congratulations on reading this far!

We have now been through a complete tour of the features needed to
implement an algebra with several realizations. The infrastructure for
realizations is not tied specifically to algebras; what we have
learned applies mutatis mutandis in full generality, for example for
implementing groups with several realizations.
