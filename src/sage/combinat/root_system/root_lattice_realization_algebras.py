r"""
Group algebras of root lattice realizations
"""
# ****************************************************************************
#       Copyright (C) 2013 Nicolas M. Thiery <nthiery at users.sf.net>
#                          Anne Schilling <anne at math.ucdavis.edu>
#                          Mark Shimozono <mshimo at vt.edu>
#                          Daniel Bump
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import functools
import operator
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.categories.algebra_functor import AlgebrasCategory
lazy_import('sage.rings.integer_ring', 'ZZ')
from sage.modules.free_module_element import vector
from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation


class Algebras(AlgebrasCategory):
    """
    The category of group algebras of root lattice realizations.

    This includes typically weight rings (group algebras of weight lattices).

    TESTS::

        sage: for ct in CartanType.samples(crystallographic=True): # long time
        ....:     TestSuite(RootSystem(ct).root_lattice().algebra(QQ)).run()
    """

    class ParentMethods:

        def _repr_(self):
            r"""
            EXAMPLES::

                sage: RootSystem(["A",2,1]).ambient_space().algebra(QQ) # indirect doctest
                Algebra of the Ambient space of the Root system of type ['A', 2, 1] over Rational Field
            """
            return "Algebra of the %s over %s" % (self.basis().keys(), self.base_ring())

        def some_elements(self):
            r"""
            Return some elements of the algebra ``self``.

            EXAMPLES::

                sage: A = RootSystem(["A",2,1]).ambient_space().algebra(QQ)
                sage: A.some_elements()
                [B[2*e[0] + 2*e[1] + 3*e[2]]...]
                sage: A.some_elements()                                                 # needs sage.graphs
                [B[2*e[0] + 2*e[1] + 3*e[2]],
                 B[-e[0] + e[2] + e['delta']],
                 B[e[0] - e[1]],
                 B[e[1] - e[2]],
                 B[e['deltacheck']],
                 B[e[0] + e['deltacheck']],
                 B[e[0] + e[1] + e['deltacheck']]]

                sage: A = RootSystem(["B",2]).weight_space().algebra(QQ)
                sage: A.some_elements()
                [B[2*Lambda[1] + 2*Lambda[2]], ... B[Lambda[1]], B[Lambda[2]]]
                sage: A.some_elements()                                                 # needs sage.graphs
                [B[2*Lambda[1] + 2*Lambda[2]],
                 B[2*Lambda[1] - 2*Lambda[2]],
                 B[-Lambda[1] + 2*Lambda[2]],
                 B[Lambda[1]],
                 B[Lambda[2]]]
            """
            return [self.monomial(weight) for weight in self.basis().keys().some_elements()]

        @cached_method
        def cartan_type(self):
            r"""
            Return the Cartan type of ``self``.

            EXAMPLES::

                sage: A = RootSystem(["A",2,1]).ambient_space().algebra(QQ)
                sage: A.cartan_type()
                ['A', 2, 1]
                sage: A = RootSystem(["B",2]).weight_space().algebra(QQ)
                sage: A.cartan_type()
                ['B', 2]
            """
            return self.basis().keys().cartan_type()

        def from_polynomial(self, p):
            """
            Construct an element of ``self`` from a polynomial `p`.

            INPUT:

            - ``p`` -- a polynomial

            EXAMPLES::

                sage: L = RootSystem(["A",2]).ambient_lattice()
                sage: KL = L.algebra(QQ)
                sage: x,y,z = QQ['x,y,z'].gens()
                sage: KL.from_polynomial(x)
                B[(1, 0, 0)]
                sage: KL.from_polynomial(x^2*y + 2*y - z)
                B[(2, 1, 0)] + 2*B[(0, 1, 0)] - B[(0, 0, 1)]

            TESTS::

                sage: KL.from_polynomial(x).leading_support().parent() is L
                True
                sage: KL.from_polynomial(x-x)
                0
                sage: KL.from_polynomial(x-x).parent() is KL
                True

            .. TODO:: make this work for Laurent polynomials too
            """
            L = self.basis().keys()
            return self.sum_of_terms((L.from_vector(vector(t)), c)
                                     for t, c in p.monomial_coefficients().items())

        @cached_method
        def divided_difference_on_basis(self, weight, i):
            r"""
            Return the result of applying the `i`-th divided difference on ``weight``.

            INPUT:

            - ``weight`` -- a weight
            - ``i`` -- an element of the index set

            .. TODO:: type free definition (Viviane's definition uses that we are in the ambient space)

            EXAMPLES::

                sage: L = RootSystem(["A",1]).ambient_space()
                sage: KL = L.algebra(QQ)
                sage: KL.divided_difference_on_basis(L((2,2)), 1) # todo: not implemented
                0
                sage: KL.divided_difference_on_basis(L((3,0)), 1) # todo: not implemented
                B[(2, 0)] + B[(1, 1)] + B[(0, 2)]
                sage: KL.divided_difference_on_basis(L((0,3)), 1) # todo: not implemented
                -B[(2, 0)] - B[(1, 1)] - B[(0, 2)]

            In type `A` and in the ambient lattice, we recover the
            usual action of divided differences polynomials::

                sage: x,y = QQ['x,y'].gens()
                sage: d = lambda p: (p - p(y,x)) / (x-y)
                sage: d(x^2*y^2)
                0
                sage: d(x^3)
                x^2 + x*y + y^2
                sage: d(y^3)
                -x^2 - x*y - y^2
            """
            raise NotImplementedError()

        @cached_method
        def isobaric_divided_difference_on_basis(self, weight, i):
            r"""
            Return the result of applying the `i`-th isobaric divided difference on ``weight``.

            INPUT:

            - ``weight`` -- a weight
            - ``i`` -- an element of the index set

            .. SEEALSO:: :meth:`demazure_operators`

            EXAMPLES::

                sage: L = RootSystem(["A",1]).ambient_space()
                sage: KL = L.algebra(QQ)
                sage: KL.isobaric_divided_difference_on_basis(L((2,2)), 1)
                B[(2, 2)]
                sage: KL.isobaric_divided_difference_on_basis(L((3,0)), 1)
                B[(1, 2)] + B[(2, 1)] + B[(3, 0)] + B[(0, 3)]
                sage: KL.isobaric_divided_difference_on_basis(L((0,3)), 1)
                -B[(1, 2)] - B[(2, 1)]

            In type `A` and in the ambient lattice, we recover the
            usual action of divided differences on polynomials::

                sage: x,y = QQ['x,y'].gens()
                sage: d = lambda p: (x*p - (x*p)(y,x)) / (x-y)
                sage: d(x^2*y^2)
                x^2*y^2
                sage: d(x^3)
                x^3 + x^2*y + x*y^2 + y^3
                sage: d(y^3)
                -x^2*y - x*y^2

            REFERENCES:

            .. [Lascoux2003] Alain Lascoux, Symmetric functions and combinatorial operators on polynomials,
               CBMS Regional Conference Series in Mathematics, 99, 2003.
            """
            P = self.basis().keys()  # the root lattice realization
            n = weight.scalar(P.simple_coroot(i))
            if n not in ZZ:
                raise ValueError("the weight does not have an integral scalar product with the coroot")
            alphai = P.simple_root(i)
            if n >= 0:
                return self.sum_of_monomials(weight-j*alphai for j in range(n + 1))
            return -self.sum_of_monomials(weight-j*alphai for j in range(n + 1, 0))

        def demazure_operators(self):
            r"""
            Return the Demazure operators acting on ``self``.

            The `i`-th Demazure operator is defined by:

            .. MATH::

                \pi_i = \frac{ 1 - e^{-\alpha_i}s_i }{ 1-e^{-\alpha_i} }

            It acts on `e^\lambda`, for `\lambda` a weight, by:

            .. MATH::

                \pi_i e^\lambda = \frac{e^\lambda - e^{-\alpha_i+s_i\lambda}}{1-e^{-\alpha_i}}

            This matches with Lascoux' definition [Lascoux2003]_ of `\pi_i`, and
            with the `i`-th Demazure operator of [Kumar1987]_, which also works for
            general Kac-Moody types.

            REFERENCES:

            .. [Kumar1987] \S. Kumar, Demazure character formula in arbitrary Kac-Moody setting,
               Invent. Math. 89 (1987), no. 2, 395-423.

            EXAMPLES:

            We compute some Schur functions, as images of dominant
            monomials under the action of the maximal isobaric divided
            difference `\Delta_{w_0}`::

                sage: L = RootSystem(["A",2]).ambient_lattice()
                sage: KL = L.algebra(QQ)
                sage: w0 = tuple(L.weyl_group().long_element().reduced_word())          # needs sage.libs.gap
                sage: pi = KL.demazure_operators()
                sage: pi0 = pi[w0]                                                      # needs sage.libs.gap
                sage: pi0(KL.monomial(L((2,1))))                                        # needs sage.libs.gap
                2*B[(1, 1, 1)] + B[(1, 2, 0)] + B[(1, 0, 2)] + B[(2, 1, 0)]
                 + B[(2, 0, 1)] + B[(0, 1, 2)] + B[(0, 2, 1)]

            Let us make the result into an actual polynomial::

                sage: P = QQ['x,y,z']
                sage: pi0(KL.monomial(L((2,1)))).expand(P.gens())                       # needs sage.libs.gap
                x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2

            This is indeed a Schur function::

                sage: s = SymmetricFunctions(QQ).s()                                    # needs sage.combinat
                sage: s[2,1].expand(3, P.variable_names())                              # needs sage.combinat
                x^2*y + x*y^2 + x^2*z + 2*x*y*z + y^2*z + x*z^2 + y*z^2

            Let us check this systematically on Schur functions of degree 6::

                sage: for p in Partitions(6, max_length=3).list():                      # needs sage.combinat sage.libs.gap
                ....:     assert (s.monomial(p).expand(3, P.variable_names())
                ....:              == pi0(KL.monomial(L(tuple(p)))).expand(P.gens()))

            We check systematically that these operators satisfy the Iwahori-Hecke algebra relations::

                sage: for cartan_type in CartanType.samples(crystallographic=True):  # long time (12s)
                ....:     L = RootSystem(cartan_type).weight_lattice()
                ....:     KL = L.algebra(QQ)
                ....:     T = KL.demazure_operators()
                ....:     T._test_relations()

                sage: L = RootSystem(['A',1,1]).weight_lattice()
                sage: KL = L.algebra(QQ)
                sage: T = KL.demazure_operators()
                sage: T._test_relations()

            .. WARNING::

                The Demazure operators are only defined if all the
                elements in the support have integral scalar products
                with the coroots (basically, they are in the weight
                lattice). Otherwise an error is raised::

                    sage: L = RootSystem(CartanType(["G",2]).dual()).ambient_space()
                    sage: KL = L.algebra(QQ)
                    sage: pi = KL.demazure_operators()
                    sage: pi[1](KL.monomial(L([0,0,1])))
                    Traceback (most recent call last):
                    ...
                    ValueError: the weight does not have an integral scalar product with the coroot
            """
            return HeckeAlgebraRepresentation(self, self.isobaric_divided_difference_on_basis, self.cartan_type(), 0, 1, side='left')

        def _test_demazure_operators(self, **options):
            """
            Test that the Demazure operators satisfy their defining formulas.

            EXAMPLES::

                sage: RootSystem(["A",2]).root_lattice().algebra(QQ)._test_demazure_operators()
            """
            tester = self._tester(**options)
            try:
                pi = self.demazure_operators()
                L = self.basis().keys()
                alpha = L.simple_roots()
                alphacheck = L.simple_coroots()
                s = L.simple_reflections()
                for i in self.cartan_type().index_set():
                    emalphai = self.monomial(-alpha[i])  # X^{-\alpha_i}
                    for weight in L.some_elements():
                        if weight.scalar(alphacheck[i]) not in ZZ:
                            # Demazure operators are not defined in this case
                            continue
                        x = self.monomial(weight)
                        result = pi[i](x)
                        tester.assertEqual(result * (self.one() - emalphai),
                                           x - emalphai * x.map_support(s[i]))
            except ImportError:
                pass

        def demazure_lusztig_operator_on_basis(self, weight, i, q1, q2, convention='antidominant'):
            r"""
            Return the result of applying the `i`-th Demazure-Lusztig operator on ``weight``.

            INPUT:

            - ``weight`` -- an element `\lambda` of the weight lattice
            - ``i`` -- an element of the index set
            - ``q1``, ``q2`` -- two elements of the ground ring
            - ``convention`` -- ``'antidominant'``, ``'bar'``, or ``'dominant'`` (default: ``'antidominant'``)

            See :meth:`demazure_lusztig_operators` for the details.

            EXAMPLES::

                sage: L = RootSystem(["A",1]).ambient_space()
                sage: K = QQ['q1,q2']
                sage: q1, q2 = K.gens()
                sage: KL = L.algebra(K)
                sage: KL.demazure_lusztig_operator_on_basis(L((2,2)), 1, q1, q2)
                q1*B[(2, 2)]
                sage: KL.demazure_lusztig_operator_on_basis(L((3,0)), 1, q1, q2)
                (q1+q2)*B[(1, 2)] + (q1+q2)*B[(2, 1)] + (q1+q2)*B[(3, 0)] + q1*B[(0, 3)]
                sage: KL.demazure_lusztig_operator_on_basis(L((0,3)), 1, q1, q2)
                (-q1-q2)*B[(1, 2)] + (-q1-q2)*B[(2, 1)] + (-q2)*B[(3, 0)]

            At `q_1=1` and `q_2=0` we recover the action of the isobaric divided differences `\pi_i`::

                sage: KL.demazure_lusztig_operator_on_basis(L((2,2)), 1, 1, 0)
                B[(2, 2)]
                sage: KL.demazure_lusztig_operator_on_basis(L((3,0)), 1, 1, 0)
                B[(1, 2)] + B[(2, 1)] + B[(3, 0)] + B[(0, 3)]
                sage: KL.demazure_lusztig_operator_on_basis(L((0,3)), 1, 1, 0)
                -B[(1, 2)] - B[(2, 1)]

            Or `1-\pi_i` for ``bar=True``::

                sage: KL.demazure_lusztig_operator_on_basis(L((2,2)), 1, 1, 0, convention='bar')
                0
                sage: KL.demazure_lusztig_operator_on_basis(L((3,0)), 1, 1, 0, convention='bar')
                -B[(1, 2)] - B[(2, 1)] - B[(0, 3)]
                sage: KL.demazure_lusztig_operator_on_basis(L((0,3)), 1, 1, 0, convention='bar')
                B[(1, 2)] + B[(2, 1)] + B[(0, 3)]

            At `q_1=1` and `q_2=-1` we recover the action of the simple reflection `s_i`::

                sage: KL.demazure_lusztig_operator_on_basis(L((2,2)), 1, 1, -1)
                B[(2, 2)]
                sage: KL.demazure_lusztig_operator_on_basis(L((3,0)), 1, 1, -1)
                B[(0, 3)]
                sage: KL.demazure_lusztig_operator_on_basis(L((0,3)), 1, 1, -1)
                B[(3, 0)]
            """
            if convention == "dominant":
                weight = -weight
            pi_on_weight = self.isobaric_divided_difference_on_basis(weight, i)
            if convention == "bar":
                pi_on_weight = self.monomial(weight) - pi_on_weight
            result = (q1+q2) * pi_on_weight - self.term(weight.simple_reflection(i), q2)
            if convention == "dominant":
                return result.map_support(operator.neg)
            else:
                return result

        def demazure_lusztig_operators(self, q1, q2, convention='antidominant'):
            r"""
            Return the Demazure-Lusztig operators acting on ``self``.

            INPUT:

            - ``q1``, ``q2`` -- two elements of the ground ring
            - ``convention`` -- "antidominant", "bar", or "dominant" (default: ``'antidominant'``)

            If `R` is the parent weight ring, the Demazure-Lusztig
            operator `T_i` is the linear map `R\rightarrow R` obtained
            by interpolating between the isobaric divided difference
            operator `\pi_i` (see :meth:`.isobaric_divided_difference_on_basis`)
            and the simple reflection `s_i`.

            .. MATH::

                (q_1+q_2) \pi_i - q_2 s_i

            The Demazure-Lusztig operators give the usual
            representation of the operator `T_i` of the (affine) Hecke
            algebra with eigenvalues `q_1` and `q_2` associated to the
            Weyl group.

            Several variants are available to match with various
            conventions used in the literature:

            - "bar" replaces `\pi_i` in the formula above by
              `\overline{\pi}_i = (1-\pi_i)`.
            - "dominant" conjugates the operator by
              `x^\lambda \mapsto x^-\lambda`.

            The names dominant and antidominant for the conventions were chosen with regards to
            the nonsymmetric Macdonald polynomials. The `Y` operators for the Macdonald polynomials
            in the "dominant" convention satisfy `Y_\lambda = T_{t_{\lambda}}` for `\lambda` dominant.
            This is also the convention used in [Haiman06]_. For the "antidominant" convention,
            `Y_\lambda = T_{t_{\lambda}}` with `\lambda` antidominant.

            .. SEEALSO::

                - :meth:`demazure_lusztig_operator_on_basis`.
                - :class:`~.non_symmetric_macdonald_polynomials.NonSymmetricMacdonaldPolynomials`.

            REFERENCES:

            .. [Lusztig1985] \G. Lusztig,
               *Equivariant K-theory and representations of Hecke algebras*,
               Proc. Amer. Math. Soc. 94 (1985), no. 2, 337-342.

            .. [Cherednik1995] \I. Cherednik,
               *Nonsymmetric Macdonald polynomials*. IMRN 10, 483-515 (1995).

            EXAMPLES::

                sage: L = RootSystem(["A",1]).ambient_space()
                sage: K = QQ['q1,q2'].fraction_field()
                sage: q1, q2 = K.gens()
                sage: KL = L.algebra(K)
                sage: T = KL.demazure_lusztig_operators(q1, q2)
                sage: Tbar = KL.demazure_lusztig_operators(q1, q2, convention='bar')
                sage: Tdominant = KL.demazure_lusztig_operators(q1, q2,
                ....:                                           convention='dominant')
                sage: x = KL.monomial(L((3,0)))
                sage: T[1](x)
                (q1+q2)*B[(1, 2)] + (q1+q2)*B[(2, 1)] + (q1+q2)*B[(3, 0)] + q1*B[(0, 3)]
                sage: Tbar[1](x)
                (-q1-q2)*B[(1, 2)] + (-q1-q2)*B[(2, 1)] + (-q1-2*q2)*B[(0, 3)]
                sage: Tbar[1](x) + T[1](x)
                (q1+q2)*B[(3, 0)] + (-2*q2)*B[(0, 3)]
                sage: Tdominant[1](x)
                (-q1-q2)*B[(1, 2)] + (-q1-q2)*B[(2, 1)] + (-q2)*B[(0, 3)]

                sage: Tdominant.Tw_inverse(1)(KL.monomial(-L.simple_root(1)))
                ((-q1-q2)/(q1*q2))*B[(0, 0)] - 1/q2*B[(1, -1)]

            We repeat similar computation in the affine setting::

                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: K = QQ['q1,q2'].fraction_field()
                sage: q1, q2 = K.gens()
                sage: KL = L.algebra(K)
                sage: T = KL.demazure_lusztig_operators(q1, q2)
                sage: Tbar = KL.demazure_lusztig_operators(q1, q2, convention='bar')
                sage: Tdominant = KL.demazure_lusztig_operators(q1, q2,
                ....:                                           convention='dominant')
                sage: e = L.basis()
                sage: x = KL.monomial(3*e[0])
                sage: T[1](x)
                (q1+q2)*B[e[0] + 2*e[1]] + (q1+q2)*B[2*e[0] + e[1]]
                 + (q1+q2)*B[3*e[0]] + q1*B[3*e[1]]
                sage: Tbar[1](x)
                (-q1-q2)*B[e[0] + 2*e[1]] + (-q1-q2)*B[2*e[0] + e[1]]
                 + (-q1-2*q2)*B[3*e[1]]
                sage: Tbar[1](x) + T[1](x)
                (q1+q2)*B[3*e[0]] + (-2*q2)*B[3*e[1]]
                sage: Tdominant[1](x)
                (-q1-q2)*B[e[0] + 2*e[1]] + (-q1-q2)*B[2*e[0] + e[1]] + (-q2)*B[3*e[1]]
                sage: Tdominant.Tw_inverse(1)(KL.monomial(-L.simple_root(1)))
                ((-q1-q2)/(q1*q2))*B[0] - 1/q2*B[e[0] - e[1]]

            One can obtain iterated operators by passing a reduced
            word or an element of the Weyl group::

                sage: T[1,2](x)
                (q1^2+2*q1*q2+q2^2)*B[e[0] + e[1] + e[2]] +
                (q1^2+2*q1*q2+q2^2)*B[e[0] + 2*e[1]] +
                (q1^2+q1*q2)*B[e[0] + 2*e[2]] + (q1^2+2*q1*q2+q2^2)*B[2*e[0] + e[1]] +
                (q1^2+q1*q2)*B[2*e[0] + e[2]] + (q1^2+q1*q2)*B[3*e[0]] +
                (q1^2+q1*q2)*B[e[1] + 2*e[2]] + (q1^2+q1*q2)*B[2*e[1] + e[2]] +
                (q1^2+q1*q2)*B[3*e[1]] + q1^2*B[3*e[2]]

            and use that to check, for example, the braid relations::

                sage: T[1,2,1](x) - T[2,1,2](x)
                0

            The operators satisfy the relations of the affine Hecke
            algebra with parameters `q_1`, `q_2`::

                sage: T._test_relations()
                sage: Tdominant._test_relations()
                sage: Tbar._test_relations() #-q2,q1+2*q2   # todo: not implemented: set the appropriate eigenvalues!

            And the `\bar{T}` are basically the inverses of the `T` s::

                sage: Tinv = KL.demazure_lusztig_operators(2/q1 + 1/q2, -1/q1,
                ....:                                      convention='bar')
                sage: [Tinv[1](T[1](x)) - x for x in KL.some_elements()]                # needs sage.graphs
                [0, 0, 0, 0, 0, 0, 0]

            We check that `\Lambda_1-\Lambda_0` is an eigenvector for
            the `Y` s in affine type::

                sage: K = QQ['q,q1,q2'].fraction_field()
                sage: q,q1,q2 = K.gens()
                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: L0 = L.classical()
                sage: Lambda = L.fundamental_weights()                                  # needs sage.graphs
                sage: alphacheck = L0.simple_coroots()
                sage: KL = L.algebra(K)
                sage: T = KL.demazure_lusztig_operators(q1, q2, convention='dominant')
                sage: Y = T.Y()
                sage: alphacheck = Y.keys().alpha()  # alpha of coroot lattice is alphacheck
                sage: alphacheck
                Finite family {0: alphacheck[0], 1: alphacheck[1], 2: alphacheck[2]}
                sage: x = KL.monomial(Lambda[1] - Lambda[0]); x                         # needs sage.graphs
                B[e[0]]

            In fact it is not exactly an eigenvector, but the extra
            '\delta` term is to be interpreted as a `q` parameter::

                sage: # needs sage.graphs
                sage: Y[alphacheck[0]](KL.one())
                q2^2/q1^2*B[0]
                sage: Y[alphacheck[1]](x)
                ((-q2^2)/(-q1^2))*B[e[0] - e['delta']]
                sage: Y[alphacheck[2]](x)
                (q1/(-q2))*B[e[0]]
                sage: KL.q_project(Y[alphacheck[1]](x),q)
                ((-q2^2)/(-q*q1^2))*B[(1, 0, 0)]

                sage: # needs sage.graphs
                sage: KL.q_project(x, q)
                B[(1, 0, 0)]
                sage: KL.q_project(Y[alphacheck[0]](x),q)
                ((-q*q1)/q2)*B[(1, 0, 0)]
                sage: KL.q_project(Y[alphacheck[1]](x),q)
                ((-q2^2)/(-q*q1^2))*B[(1, 0, 0)]
                sage: KL.q_project(Y[alphacheck[2]](x),q)
                (q1/(-q2))*B[(1, 0, 0)]

            We now check systematically that the Demazure-Lusztig
            operators satisfy the relations of the Iwahori-Hecke
            algebra::

                sage: K = QQ['q1,q2']
                sage: q1, q2 = K.gens()
                sage: for cartan_type in CartanType.samples(crystallographic=True):  # long time (12s)
                ....:    L = RootSystem(cartan_type).root_lattice()
                ....:    KL = L.algebra(K)
                ....:    T = KL.demazure_lusztig_operators(q1,q2)
                ....:    T._test_relations()

                sage: for cartan_type in CartanType.samples(crystallographic=True):  # long time (12s)
                ....:    L = RootSystem(cartan_type).weight_lattice()
                ....:    KL = L.algebra(K)
                ....:    T = KL.demazure_lusztig_operators(q1,q2)
                ....:    T._test_relations()

            Recall that the Demazure-Lusztig operators are only
            defined when all monomials belong to the weight lattice.
            Thus, in the group algebra of the ambient space, we need
            to specify explicitly the elements on which to run the
            tests::

                sage: for cartan_type in CartanType.samples(crystallographic=True):  # long time (12s)
                ....:    L = RootSystem(cartan_type).ambient_space()
                ....:    KL = L.algebra(K)
                ....:    weight_lattice = RootSystem(cartan_type).weight_lattice(extended=L.is_extended())
                ....:    elements = [ KL.monomial(L(weight)) for weight in weight_lattice.some_elements() ]
                ....:    T = KL.demazure_lusztig_operators(q1,q2)
                ....:    T._test_relations(elements=elements)
            """
            T_on_basis = functools.partial(self.demazure_lusztig_operator_on_basis,
                                           q1=q1, q2=q2, convention=convention)
            return HeckeAlgebraRepresentation(self, T_on_basis, self.cartan_type(), q1, q2, side='left')

        def demazure_lusztig_operator_on_classical_on_basis(self, weight, i, q, q1, q2, convention='antidominant'):
            r"""
            Return the result of applying the `i`-th Demazure-Lusztig operator on the classical weight ``weight`` embedded at level 0.

            INPUT:

            - ``weight`` -- a classical weight `\lambda`
            - ``i`` -- an element of the index set
            - ``q1``, ``q2`` -- two elements of the ground ring
            - ``convention`` -- ``'antidominant'``, ``'bar'``, or ``'dominant'`` (default: ``'antidominant'``)

            See :meth:`demazure_lusztig_operators` for the details.

            .. TODO::

                - Optimize the code to only do the embedding/projection for T_0
                - Add an option to specify at which level one wants to
                  work. Currently this is level 0.

            EXAMPLES::

                sage: L = RootSystem(["A",1,1]).ambient_space()
                sage: L0 = L.classical()
                sage: K = QQ['q,q1,q2']
                sage: q, q1, q2 = K.gens()
                sage: KL = L.algebra(K)
                sage: KL0 = L0.algebra(K)

            These operators coincide with the usual Demazure-Lusztig
            operators::

                sage: KL.demazure_lusztig_operator_on_classical_on_basis(L0((2,2)),     # needs sage.graphs
                ....:                                                    1, q, q1, q2)
                q1*B[(2, 2)]
                sage: KL0.demazure_lusztig_operator_on_basis(L0((2,2)), 1, q1, q2)
                q1*B[(2, 2)]

                sage: KL.demazure_lusztig_operator_on_classical_on_basis(L0((3,0)),     # needs sage.graphs
                ....:                                                    1, q, q1, q2)
                (q1+q2)*B[(1, 2)] + (q1+q2)*B[(2, 1)] + (q1+q2)*B[(3, 0)] + q1*B[(0, 3)]
                sage: KL0.demazure_lusztig_operator_on_basis(L0((3,0)), 1, q1, q2)
                (q1+q2)*B[(1, 2)] + (q1+q2)*B[(2, 1)] + (q1+q2)*B[(3, 0)] + q1*B[(0, 3)]

            except that we now have an action of `T_0`, which introduces some `q` s::

                sage: KL.demazure_lusztig_operator_on_classical_on_basis(L0((2,2)),     # needs sage.graphs
                ....:                                                    0, q, q1, q2)
                q1*B[(2, 2)]
                sage: KL.demazure_lusztig_operator_on_classical_on_basis(L0((3,0)),     # needs sage.graphs
                ....:                                                    0, q, q1, q2)
                (-q^2*q1-q^2*q2)*B[(1, 2)] + (-q*q1-q*q2)*B[(2, 1)] + (-q^3*q2)*B[(0, 3)]
            """
            L = self.basis().keys()
            weight = L.embed_at_level(weight, 0)
            return self.q_project(self.demazure_lusztig_operator_on_basis(weight, i, q1, q2, convention=convention), q)

        def demazure_lusztig_operators_on_classical(self, q, q1, q2, convention='antidominant'):
            r"""
            Return the Demazure-Lusztig operators acting at level 1 on ``self.classical()``.

            INPUT:

            - ``q``, ``q1``, ``q2`` -- three elements of the ground ring
            - ``convention`` -- ``'antidominant'``, ``'bar'``, or ``'dominant'`` (default: ``'antidominant'``)

            Let `KL` be the group algebra of an affine weight lattice
            realization `L`. The Demazure-Lusztig operators for `KL`
            act on the group algebra of the corresponding classical
            weight lattice by embedding it at level 1, and projecting
            back.

            .. SEEALSO::

                - :meth:`demazure_lusztig_operators`.
                - :meth:`demazure_lusztig_operator_on_classical_on_basis`.
                - :meth:`q_project`

            EXAMPLES::

                sage: L = RootSystem(["A",1,1]).ambient_space()
                sage: K = QQ['q,q1,q2'].fraction_field()
                sage: q, q1, q2 = K.gens()
                sage: KL = L.algebra(K)
                sage: KL0 = KL.classical()
                sage: L0 = KL0.basis().keys()
                sage: T = KL.demazure_lusztig_operators_on_classical(q, q1, q2)         # needs sage.graphs

                sage: x = KL0.monomial(L0((3,0))); x
                B[(3, 0)]

            For `T_1,\dots` we recover the usual Demazure-Lusztig operators::

                sage: T[1](x)                                                           # needs sage.graphs
                (q1+q2)*B[(1, 2)] + (q1+q2)*B[(2, 1)] + (q1+q2)*B[(3, 0)] + q1*B[(0, 3)]

            For `T_0`, we can note that, in the projection, `\delta`
            is mapped to `q`::

                sage: T[0](x)                                                           # needs sage.graphs
                (-q^2*q1-q^2*q2)*B[(1, 2)] + (-q*q1-q*q2)*B[(2, 1)] + (-q^3*q2)*B[(0, 3)]

            Note that there is no translation part, and in particular
            1 is an eigenvector for all `T_i`'s::

                sage: # needs sage.graphs
                sage: T[0](KL0.one())
                q1*B[(0, 0)]
                sage: T[1](KL0.one())
                q1*B[(0, 0)]
                sage: Y = T.Y()
                sage: alphacheck = Y.keys().simple_roots()
                sage: Y[alphacheck[0]](KL0.one())
                ((-q2)/(q*q1))*B[(0, 0)]

            Matching with Ion Bogdan's hand calculations from 3/15/2013::

                sage: L = RootSystem(["A",1,1]).weight_space(extended=True)
                sage: K = QQ['q,u'].fraction_field()
                sage: q, u = K.gens()
                sage: KL = L.algebra(K)
                sage: KL0 = KL.classical()
                sage: L0 = KL0.basis().keys()
                sage: omega = L0.fundamental_weights()

                sage: # needs sage.graphs
                sage: T = KL.demazure_lusztig_operators_on_classical(q, u, -1/u,
                ....:                                                convention='dominant')
                sage: Y = T.Y()
                sage: alphacheck = Y.keys().simple_roots()
                sage: Ydelta = Y[Y.keys().null_root()]
                sage: Ydelta.word, Ydelta.signs, Ydelta.scalar
                ((), (), 1/q)
                sage: Y1 = Y[alphacheck[1]]
                sage: Y1.word, Y1.signs, Y1.scalar  # This is T_0 T_1 (T_1 acts first, then T_0); Ion gets T_1 T_0
                ((1, 0), (1, 1), 1)
                sage: Y0 = Y[alphacheck[0]]
                sage: Y0.word, Y0.signs, Y0.scalar  # This is 1/q T_1^-1 T_0^-1
                ((0, 1), (-1, -1), 1/q)

            Note that the following computations use the "dominant" convention::

                sage: # needs sage.graphs
                sage: T0 = T.Tw(0)
                sage: T0(KL0.monomial(omega[1]))
                q*u*B[-Lambda[1]] + ((u^2-1)/u)*B[Lambda[1]]
                sage: T0(KL0.monomial(2*omega[1]))
                ((q*u^2-q)/u)*B[0] + q^2*u*B[-2*Lambda[1]] + ((u^2-1)/u)*B[2*Lambda[1]]
                sage: T0(KL0.monomial(-omega[1]))
                1/(q*u)*B[Lambda[1]]
                sage: T0(KL0.monomial(-2*omega[1]))
                ((-u^2+1)/(q*u))*B[0] + 1/(q^2*u)*B[2*Lambda[1]]
            """
            # In type BC dual we used q^2 and q elsewhere
            # Not sure this is the right thing to do or just a workaround ...
            # This probably makes up for the fact that, in type BC
            # dual, the null coroot is twice Sage's deltacheck
            # whereas the null root is delta. So we need to map delta
            # to q^2 in the q_projection.
            # Should this go in q_project instead?
            ct = self.cartan_type()
            a0check = ct.acheck()[ct.special_node()]
            T_on_basis = functools.partial(self.demazure_lusztig_operator_on_classical_on_basis,
                                           q1=q1, q2=q2, q=q**a0check, convention=convention)
            return HeckeAlgebraRepresentation(self.classical(), T_on_basis, self.cartan_type(), q1=q1, q2=q2, q=q, side='left')

        @cached_method
        def T0_check_on_basis(self, q1, q2, convention='antidominant'):
            r"""
            Return the `T_0^\vee` operator acting on the basis.

            This implements the formula for `T_{0'}` in Section 6.12 of [Haiman06]_.

            REFERENCES:

            .. [Haiman06] \M. Haiman, Cherednik algebras, Macdonald polynomials and combinatorics, ICM 2006.

            .. WARNING::

                The current implementation probably returns just
                nonsense, if the convention is not "dominant".

            EXAMPLES::

                sage: K = QQ['q1,q2'].fraction_field()
                sage: q1,q2 = K.gens()

                sage: L = RootSystem(["A",1,1]).ambient_space()
                sage: L0 = L.classical()
                sage: KL = L.algebra(K)
                sage: some_weights = L.fundamental_weights()                            # needs sage.graphs
                sage: f = KL.T0_check_on_basis(q1,q2, convention='dominant')            # needs sage.graphs
                sage: f(L0.zero())                                                      # needs sage.graphs
                (q1+q2)*B[(0, 0)] + q1*B[(1, -1)]

                sage: L = RootSystem(["A",3,1]).ambient_space()
                sage: L0 = L.classical()
                sage: KL = L.algebra(K)
                sage: some_weights = L0.fundamental_weights()
                sage: f = KL.T0_check_on_basis(q1,q2, convention='dominant')            # needs sage.graphs
                sage: f(L0.zero())       # not checked                                  # needs sage.graphs
                (q1+q2)*B[(0, 0, 0, 0)] + q1^3/q2^2*B[(1, 0, 0, -1)]

            The following results have not been checked::

                sage: for x in some_weights:                                            # needs sage.graphs
                ....:     print("{} : {}".format(x, f(x)))
                (1, 0, 0, 0) : q1*B[(1, 0, 0, 0)]
                (1, 1, 0, 0) : q1*B[(1, 1, 0, 0)]
                (1, 1, 1, 0) : q1*B[(1, 1, 1, 0)]

            Some examples for type `B_2^{(1)}` dual::

                sage: L = RootSystem("B2~*").ambient_space()
                sage: L0 = L.classical()
                sage: e = L.basis()
                sage: K = QQ['q,u'].fraction_field()
                sage: q,u = K.gens()
                sage: q1 = u
                sage: q2 = -1/u
                sage: KL = L.algebra(K)
                sage: KL0 = KL.classical()
                sage: f = KL.T0_check_on_basis(q1,q2, convention='dominant')            # needs sage.graphs
                sage: T = KL.twisted_demazure_lusztig_operators(q1,q2, convention='dominant')

            Direct calculation::

                sage: T.Tw(0)(KL0.monomial(L0([0,0])))                                  # needs sage.graphs
                ((u^2-1)/u)*B[(0, 0)] + u^3*B[(1, 1)]
                sage: KL.T0_check_on_basis(q1,q2, convention='dominant')(L0([0,0]))     # needs sage.graphs
                ((u^2-1)/u)*B[(0, 0)] + u^3*B[(1, 1)]

            Step by step calculation, comparing by hand with Mark Shimozono::

                sage: res = T.Tw(2)(KL0.monomial(L0([0,0]))); res
                u*B[(0, 0)]
                sage: res = res * KL0.monomial(L0([-1,1])); res
                u*B[(-1, 1)]
                sage: res = T.Tw_inverse(1)(res); res
                (u^2-1)*B[(0, 0)] + u^2*B[(1, -1)]
                sage: res = T.Tw_inverse(2)(res); res
                ((u^2-1)/u)*B[(0, 0)] + u^3*B[(1, 1)]
            """
            L = self.basis().keys()
            ct = L.cartan_type()
            special_node = ct.special_node()
            a0 = ct.a()[special_node]
            A0 = self.classical()
            T = A0.demazure_lusztig_operators(q1, q2, convention=convention)
            # TODO: use the formula expressing the inverse of T as a Demazure Lusztig operator? Or go through the affine action of T_0 for the dual
            L0 = A0.basis().keys()
            # The dominant short root of the classical system
            if ct.type() == 'BC':
                # CHECKME: this is not exactly phi, but phi rescaled
                # appropriately so that it's in the orbit of the
                # simple classical roots
                phi = -a0*L0(L.simple_roots()[0])
            else:
                phi = L0(L0.root_system.coroot_lattice().highest_root().associated_coroot())
            # Variant: try to fetch it from the other affinization; something like:
            # The a0 only has an influence in type BC; it handles the fact that alpha_0
            # is not in the orbit of the classical roots
            #phi1 = - L0(L'.other_affinization().simple_roots()[special_node]) * a0
            #assert phi == phi1

            j, v = phi.to_simple_root(reduced_word=True)
            translation = A0.monomial(-L0.simple_root(j)/a0)
            Tv = T[v]
            Tinv = T.Tw_inverse(v+(j,))

            def T0_check(weight):
                return -q1*q2*Tinv( translation * Tv(A0.monomial(weight)))
            # For debugging purposes
            T0_check.phi = phi
            T0_check.j = j
            T0_check.v = v
            return T0_check

        @cached_method
        def classical(self):
            """
            Return the group algebra of the corresponding classical lattice.

            EXAMPLES::

                sage: KL = RootSystem(["A",2,1]).ambient_space().algebra(QQ)
                sage: KL.classical()
                Algebra of the Ambient space of the Root system of type ['A', 2] over Rational Field
            """
            return self.basis().keys().classical().algebra(self.base_ring())

        def q_project_on_basis(self, l, q):
            r"""
            Return the monomial `c * cl(l)`  in the group algebra of the classical lattice.

            INPUT:

            - ``l`` -- an element of the root lattice realization
            - ``q`` -- an element of the ground ring

            Here, `cl(l)` is the projection of `l` in the classical
            lattice, and `c` is the coefficient of `l` in `\delta`.

            .. SEEALSO:: :meth:`q_project_on_basis`

            EXAMPLES::

                sage: K = QQ['q'].fraction_field()
                sage: q = K.gen()
                sage: KL = RootSystem(["A",2,1]).ambient_space().algebra(K)
                sage: L = KL.basis().keys()
                sage: e = L.basis()
                sage: KL.q_project_on_basis(4*e[1] + 3*e[2]
                ....:                        + e['deltacheck'] - 2*e['delta'], q)
                1/q^2*B[(0, 4, 3)]
            """
            KL0 = self.classical()
            L0 = KL0.basis().keys()
            return KL0.term(L0(l), q**l["delta"])

        def q_project(self, x, q):
            r"""
            Implement the `q`-projection morphism from ``self`` to the group algebra of the classical space.

            INPUT:

            - ``x`` -- an element of the group algebra of ``self``
            - ``q`` -- an element of the ground ring

            This is an algebra morphism mapping `\delta` to `q` and
            `X^b` to its classical counterpart for the other elements
            `b` of the basis of the realization.

            EXAMPLES::

                sage: K = QQ['q'].fraction_field()
                sage: q = K.gen()
                sage: KL = RootSystem(["A",2,1]).ambient_space().algebra(K)
                sage: L = KL.basis().keys()
                sage: e = L.basis()
                sage: x = (KL.an_element()
                ....:       + KL.monomial(4*e[1] + 3*e[2]
                ....:                      + e['deltacheck'] - 2*e['delta'])); x
                B[2*e[0] + 2*e[1] + 3*e[2]]
                 + B[4*e[1] + 3*e[2] - 2*e['delta'] + e['deltacheck']]
                sage: KL.q_project(x, q)
                B[(2, 2, 3)] + 1/q^2*B[(0, 4, 3)]

                sage: KL = RootSystem(["BC",3,2]).ambient_space().algebra(K)
                sage: L = KL.basis().keys()
                sage: e = L.basis()
                sage: x = (KL.an_element()
                ....:       + KL.monomial(4*e[1] + 3*e[2]
                ....:                      + e['deltacheck'] - 2*e['delta'])); x
                B[2*e[0] + 2*e[1] + 3*e[2]]
                 + B[4*e[1] + 3*e[2] - 2*e['delta'] + e['deltacheck']]
                sage: KL.q_project(x, q)
                B[(2, 2, 3)] + 1/q^2*B[(0, 4, 3)]

            .. WARNING::

                Recall that the null root, usually denoted `\delta`,
                is in fact ``a[0]\delta`` in Sage's notation, in order
                to avoid half integer coefficients (this only makes a
                difference in type BC). Similarly, what's usually
                denoted `q` is in fact ``q^a[0]`` in Sage's notations,
                to avoid manipulating square roots::

                    sage: KL.q_project(KL.monomial(L.null_root()),q)                    # needs sage.graphs
                    q^2*B[(0, 0, 0)]
            """
            L0 = self.classical()
            return L0.linear_combination( (self.q_project_on_basis(l, q), c) for l,c in x )

        def twisted_demazure_lusztig_operator_on_basis(self, weight, i, q1, q2, convention='antidominant'):
            r"""
            Return the twisted Demazure-Lusztig operator acting on the basis.

            INPUT:

            - ``weight`` -- an element `\lambda` of the weight lattice
            - ``i`` -- an element of the index set
            - ``q1``, ``q2`` -- two elements of the ground ring
            - ``convention`` -- ``'antidominant'``, ``'bar'``, or ``'dominant'`` (default: ``'antidominant'``)

            .. SEEALSO:: :meth:`twisted_demazure_lusztig_operators`

            EXAMPLES::

                sage: L = RootSystem(["A",3,1]).ambient_space()
                sage: e = L.basis()
                sage: K = QQ['q1,q2'].fraction_field()
                sage: q1, q2 = K.gens()
                sage: KL = L.algebra(K)
                sage: Lambda = L.classical().fundamental_weights()
                sage: KL.twisted_demazure_lusztig_operator_on_basis(
                ....:     Lambda[1] + 2*Lambda[2], 1, q1, q2, convention='dominant')
                (-q2)*B[(2, 3, 0, 0)]
                sage: KL.twisted_demazure_lusztig_operator_on_basis(
                ....:     Lambda[1] + 2*Lambda[2], 2, q1, q2, convention='dominant')
                (-q1-q2)*B[(3, 1, 1, 0)] + (-q2)*B[(3, 0, 2, 0)]
                sage: KL.twisted_demazure_lusztig_operator_on_basis(
                ....:     Lambda[1] + 2*Lambda[2], 3, q1, q2, convention='dominant')
                q1*B[(3, 2, 0, 0)]
                sage: KL.twisted_demazure_lusztig_operator_on_basis(                    # needs sage.graphs
                ....:     Lambda[1]+2*Lambda[2], 0, q1, q2, convention='dominant')
                ((q1*q2+q2^2)/q1)*B[(1, 2, 1, 1)] + ((q1*q2+q2^2)/q1)*B[(1, 2, 2, 0)]
                 + q2^2/q1*B[(1, 2, 0, 2)] + ((q1^2+2*q1*q2+q2^2)/q1)*B[(2, 1, 1, 1)]
                 + ((q1^2+2*q1*q2+q2^2)/q1)*B[(2, 1, 2, 0)]
                 + ((q1*q2+q2^2)/q1)*B[(2, 1, 0, 2)]
                 + ((q1^2+2*q1*q2+q2^2)/q1)*B[(2, 2, 1, 0)]
                 + ((q1*q2+q2^2)/q1)*B[(2, 2, 0, 1)]
            """
            if i == 0: # should use the special node
                if convention != "dominant":
                    raise NotImplementedError("The twisted Demazure-Lusztig operator T_0 is only implemented in the dominant convention")
                return self.T0_check_on_basis(q1, q2, convention=convention)(weight)
            else:
                L = self.classical()
                return L.demazure_lusztig_operators(q1, q2, convention=convention)[i](L.monomial(weight))

        def twisted_demazure_lusztig_operators(self, q1, q2, convention='antidominant'):
            r"""
            Return the twisted Demazure-Lusztig operators acting on ``self``.

            INPUT:

            - ``q1``, ``q2`` -- two elements of the ground ring
            - ``convention`` -- ``'antidominant'``, ``'bar'``, or ``'dominant'`` (default: ``'antidominant'``)

            .. WARNING::

                - the code is currently only tested for `q_1q_2=-1`
                - only the ``'dominant'`` convention is functional for `i=0`

            For `T_1,\ldots,T_n`, these operators are the usual
            Demazure-Lusztig operators. On the other hand, the
            operator `T_0` is twisted::

                sage: L = RootSystem(["A",3,1]).ambient_space()
                sage: e = L.basis()
                sage: K = QQ['q1,q2'].fraction_field()
                sage: q1, q2 = K.gens()
                sage: KL = L.algebra(K)
                sage: T = KL.twisted_demazure_lusztig_operators(q1, q2, convention='dominant')
                sage: T._test_relations()

            TESTS:

            The following computations were checked with Mark Shimozono for type `A_1^{(1)}`::

                sage: L = RootSystem(["A",1,1]).ambient_space()
                sage: e = L.basis()
                sage: K = QQ['q1,q2'].fraction_field()
                sage: q1,q2 = K.gens()
                sage: KL = L.algebra(K)
                sage: T = KL.twisted_demazure_lusztig_operators(q1, q2, convention='dominant')
                sage: T._test_relations()
                sage: L0 = L.classical()
                sage: alpha = L0.simple_roots()
                sage: T.Ti_on_basis(L0.zero(), 1)
                q1*B[(0, 0)]
                sage: T.Ti_inverse_on_basis(L0.zero(), 1)
                1/q1*B[(0, 0)]
                sage: T.Ti_on_basis(alpha[1], 1)
                (-q1-q2)*B[(0, 0)] + (-q2)*B[(-1, 1)]
                sage: T.Ti_inverse_on_basis(alpha[1], 1)
                ((q1+q2)/(q1*q2))*B[(0, 0)] + 1/q1*B[(-1, 1)] + ((q1+q2)/(q1*q2))*B[(1, -1)]
                sage: T.Ti_on_basis(L0.zero(), 0)                                       # needs sage.graphs
                (q1+q2)*B[(0, 0)] + q1*B[(1, -1)]

            The next computations were checked with Mark Shimozono for type `A_2^{(1)}`::

                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: e = L.basis()
                sage: K = QQ['u'].fraction_field()
                sage: u = K.gen()
                sage: KL = L.algebra(K)
                sage: T = KL.twisted_demazure_lusztig_operators(u, -~u, convention='dominant')
                sage: T._test_relations()
                sage: L0 = L.classical()
                sage: KL0 = L0.algebra(K)
                sage: alpha = L0.simple_roots()

                sage: phi = L0.highest_root(); phi
                (1, 0, -1)
                sage: phi.to_simple_root(reduced_word=True)
                (2, (1,))
                sage: res = T.Ti_on_basis(L0([1,0,1]), 1); res
                1/u*B[(0, 1, 1)]
                sage: res = res * KL0.monomial(-alpha[2]); res
                1/u*B[(0, 0, 2)]
                sage: res = T.Tw_inverse(2)(res); res
                ((u^2-1)/u^2)*B[(0, 1, 1)] + B[(0, 2, 0)]
                sage: res = T.Tw_inverse(1)(res); res
                ((u^2-1)/u)*B[(1, 1, 0)] + ((u^2-1)/u)*B[(1, 0, 1)] + u*B[(2, 0, 0)]

            .. TODO::

                Choose a good set of Cartan Type to run on. Rank >4 is
                too big. But `C_1` and `B_1` are boring.

            We now check systematically that those operators satisfy
            the relations of the Iwahori-Hecke algebra::

                sage: K = QQ['q1,q2'].fraction_field()
                sage: q1, q2 = K.gens()
                sage: for cartan_type in CartanType.samples(affine=True, crystallographic=True):  # long time (12s)
                ....:     if cartan_type.rank() > 4: continue
                ....:     if cartan_type.type() == 'BC': continue
                ....:     KL = RootSystem(cartan_type).weight_lattice().algebra(K)
                ....:     T = KL.twisted_demazure_lusztig_operators(q1, q2, convention='dominant')
                ....:     T._test_relations()

            .. TODO::

                Investigate why `T_0^\vee` currently does not satisfy
                the quadratic relation in type `BC`. This should
                hopefully be fixed when `T_0^\vee` will have a more
                uniform implementation::

                    sage: cartan_type = CartanType(["BC",1,2])
                    sage: KL = RootSystem(cartan_type).weight_lattice().algebra(K)
                    sage: T = KL.twisted_demazure_lusztig_operators(q1,q2, convention='dominant')
                    sage: T._test_relations()                                           # needs sage.graphs
                    Traceback (most recent call last):
                    ... tester.assertTrue(Ti(Ti(x,i,-q2),i,-q1).is_zero()) ...
                    AssertionError: False is not true

            Comparison with T0::

                sage: L = RootSystem(["A",2,1]).ambient_space()
                sage: e = L.basis()
                sage: K = QQ['t,q'].fraction_field()
                sage: t,q = K.gens()
                sage: q1 = t
                sage: q2 = -1
                sage: KL = L.algebra(K)
                sage: L0 = L.classical()
                sage: T = KL.demazure_lusztig_operators(q1,q2, convention='dominant')
                sage: def T0(*l0): return KL.q_project(T[0].on_basis()(L.embed_at_level(L0(l0), 1)), q)
                sage: T0_check_on_basis = KL.T0_check_on_basis(q1, q2,                  # needs sage.graphs
                ....:                                          convention='dominant')
                sage: def T0c(*l0): return T0_check_on_basis(L0(l0))

                sage: T0(0,0,1)                                 # not double checked    # needs sage.graphs
                ((-t+1)/q)*B[(1, 0, 0)] + 1/q^2*B[(2, 0, -1)]
                sage: T0c(0,0,1)                                                        # needs sage.graphs
                (t^2-t)*B[(1, 0, 0)] + (t^2-t)*B[(1, 1, -1)] + t^2*B[(2, 0, -1)] + (t-1)*B[(0, 0, 1)]
            """
            T_on_basis = functools.partial(self.twisted_demazure_lusztig_operator_on_basis,
                                           q1=q1, q2=q2, convention=convention)
            return HeckeAlgebraRepresentation(self.classical(),
                                              T_on_basis,
                                              self.cartan_type().classical().dual().affine().dual(),
                                              q1, q2,
                                              side='left')

    class ElementMethods:

        def acted_upon(self, w):
            """
            Implement the action of ``w`` on ``self``.

            INPUT:

            - ``w`` -- an element of the Weyl group acting on the underlying weight lattice realization

            EXAMPLES::

                sage: L = RootSystem(["A",3]).ambient_space()
                sage: W = L.weyl_group()                                                # needs sage.libs.gap
                sage: M = L.algebra(QQ['q','t'])
                sage: m = M.an_element(); m  # TODO: investigate why we don't get something more interesting
                B[(2, 2, 3, 0)]
                sage: m = (m+1)^2; m
                B[(0, 0, 0, 0)] + 2*B[(2, 2, 3, 0)] + B[(4, 4, 6, 0)]
                sage: w = W.an_element(); w.reduced_word()                              # needs sage.libs.gap
                [1, 2, 3]
                sage: m.acted_upon(w)                                                   # needs sage.libs.gap
                B[(0, 0, 0, 0)] + 2*B[(0, 2, 2, 3)] + B[(0, 4, 4, 6)]
            """
            return self.map_support(w.action)

        def expand(self, alphabet):
            """
            Expand ``self`` into variables in the ``alphabet``.

            INPUT:

            - ``alphabet`` -- a non empty list/tuple of (invertible) variables in a ring to expand in

            EXAMPLES::

                sage: L = RootSystem(["A",2]).ambient_lattice()
                sage: KL = L.algebra(QQ)
                sage: p = KL.an_element() + KL.sum_of_monomials(L.some_elements()); p
                B[(1, 0, 0)] + B[(1, -1, 0)] + B[(1, 1, 0)] + 2*B[(2, 2, 3)] + B[(0, 1, -1)]
                sage: F = LaurentPolynomialRing(QQ, 'x,y,z')
                sage: p.expand(F.gens())
                2*x^2*y^2*z^3 + x*y + x + y*z^-1 + x*y^-1

            TESTS::

                sage: type(p.expand(F.gens()))
                <class 'sage.rings.polynomial.laurent_polynomial_mpair.LaurentPolynomial_mpair'>

                sage: p = KL.zero()
                sage: p.expand(F.gens())
                0
                sage: type(p.expand(F.gens()))
                <class 'sage.rings.polynomial.laurent_polynomial_mpair.LaurentPolynomial_mpair'>
            """
            codomain = alphabet[0].parent()
            return codomain.sum(c * prod(X**int(n)
                                         for X, n in zip(alphabet, vector(m)))
                                for m, c in self)
