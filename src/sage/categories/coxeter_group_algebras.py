r"""
Coxeter Group Algebras
"""
import functools
from sage.misc.cachefunc import cached_method
from sage.categories.algebra_functor import AlgebrasCategory

class CoxeterGroupAlgebras(AlgebrasCategory):

    class ParentMethods:

        def demazure_lusztig_operator_on_basis(self, w, i, q1, q2,
                                               side="right"):
            r"""
            Return the result of applying the `i`-th Demazure Lusztig
            operator on ``w``.

            INPUT:

            - ``w`` -- an element of the Coxeter group
            - ``i`` -- an element of the index set
            - ``q1,q2`` -- two elements of the ground ring
            - ``bar`` -- a boolean (default ``False``)

            See :meth:`demazure_lusztig_operators` for details.

            EXAMPLES::

                sage: W = WeylGroup(["B",3])
                sage: W.element_class._repr_=lambda x: "".join(str(i) for i in x.reduced_word())
                sage: K = QQ['q1,q2']
                sage: q1, q2 = K.gens()
                sage: KW = W.algebra(K)
                sage: w = W.an_element()
                sage: KW.demazure_lusztig_operator_on_basis(w, 0, q1, q2)
                (-q2)*323123 + (q1+q2)*123
                sage: KW.demazure_lusztig_operator_on_basis(w, 1, q1, q2)
                q1*1231
                sage: KW.demazure_lusztig_operator_on_basis(w, 2, q1, q2)
                q1*1232
                sage: KW.demazure_lusztig_operator_on_basis(w, 3, q1, q2)
                (q1+q2)*123 + (-q2)*12

            At `q_1=1` and `q_2=0` we recover the action of the
            isobaric divided differences `\pi_i`::

                sage: KW.demazure_lusztig_operator_on_basis(w, 0, 1, 0)
                123
                sage: KW.demazure_lusztig_operator_on_basis(w, 1, 1, 0)
                1231
                sage: KW.demazure_lusztig_operator_on_basis(w, 2, 1, 0)
                1232
                sage: KW.demazure_lusztig_operator_on_basis(w, 3, 1, 0)
                123

            At `q_1=1` and `q_2=-1` we recover the action of the
            simple reflection `s_i`::

                sage: KW.demazure_lusztig_operator_on_basis(w, 0, 1, -1)
                323123
                sage: KW.demazure_lusztig_operator_on_basis(w, 1, 1, -1)
                1231
                sage: KW.demazure_lusztig_operator_on_basis(w, 2, 1, -1)
                1232
                sage: KW.demazure_lusztig_operator_on_basis(w, 3, 1, -1)
                12
            """
            return (q1+q2) * self.monomial(w.apply_simple_projection(i,side=side)) - self.term(w.apply_simple_reflection(i, side=side), q2)

        def demazure_lusztig_operators(self, q1, q2, side="right", affine=True):
            r"""
            Return the Demazure Lusztig operators acting on ``self``.

            INPUT:

            - ``q1,q2`` -- two elements of the ground ring `K`
            - ``side`` -- ``"left"`` or ``"right"`` (default: ``"right"``);
              which side to act upon
            - ``affine`` -- a boolean (default: ``True``)

            The Demazure-Lusztig operator `T_i` is the linear map
            `R \to R` obtained by interpolating between the
            simple projection `\pi_i` (see
            :meth:`CoxeterGroups.ElementMethods.simple_projection`)
            and the simple reflection `s_i` so that `T_i` has
            eigenvalues `q_1` and `q_2`:

            .. MATH::

                (q_1 + q_2) \pi_i - q_2 s_i.

            The Demazure-Lusztig operators give the usual
            representation of the operators `T_i` of the `q_1,q_2`
            Hecke algebra associated to the Coxeter group.

            For a finite Coxeter group, and if ``affine=True``, the
            Demazure-Lusztig operators `T_1,\dots,T_n` are completed
            by `T_0` to implement the level `0` action of the affine
            Hecke algebra.

            EXAMPLES::

                sage: W = WeylGroup(["B",3])
                sage: W.element_class._repr_=lambda x: "".join(str(i) for i in x.reduced_word())
                sage: K = QQ['q1,q2']
                sage: q1, q2 = K.gens()
                sage: KW = W.algebra(K)
                sage: T = KW.demazure_lusztig_operators(q1, q2, affine=True)
                sage: x = KW.monomial(W.an_element()); x
                123
                sage: T[0](x)
                (-q2)*323123 + (q1+q2)*123
                sage: T[1](x)
                q1*1231
                sage: T[2](x)
                q1*1232
                sage: T[3](x)
                (q1+q2)*123 + (-q2)*12

                sage: T._test_relations()

            .. NOTE::

                For a finite Weyl group `W`, the level 0 action of the
                affine Weyl group `\tilde W` only depends on the
                Coxeter diagram of the affinization, not its Dynkin
                diagram. Hence it is possible to explore all cases
                using only untwisted affinizations.
            """
            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            W = self.basis().keys()
            cartan_type = W.cartan_type()
            if affine and cartan_type.is_finite():
                cartan_type = cartan_type.affine()
            T_on_basis = functools.partial(self.demazure_lusztig_operator_on_basis, q1=q1, q2=q2, side=side)
            return HeckeAlgebraRepresentation(self, T_on_basis, cartan_type, q1, q2)

        @cached_method
        def demazure_lusztig_eigenvectors(self, q1, q2):
            r"""
            Return the family of eigenvectors for the Cherednik operators.

            INPUT:

            - ``self`` -- a finite Coxeter group `W`
            - ``q1,q2`` -- two elements of the ground ring `K`

            The affine Hecke algebra `H_{q_1,q_2}(\tilde W)` acts on
            the group algebra of `W` through the Demazure-Lusztig
            operators `T_i`. Its Cherednik operators `Y^\lambda` can
            be simultaneously diagonalized as long as `q_1/q_2` is not
            a small root of unity [HST2008]_.

            This method returns the family of joint eigenvectors,
            indexed by `W`.

            .. SEEALSO::

                - :meth:`demazure_lusztig_operators`
                - :class:`sage.combinat.root_system.hecke_algebra_representation.CherednikOperatorsEigenvectors`

            EXAMPLES::

                sage: W = WeylGroup(["B",2])
                sage: W.element_class._repr_=lambda x: "".join(str(i) for i in x.reduced_word())
                sage: K = QQ['q1,q2'].fraction_field()
                sage: q1, q2 = K.gens()
                sage: KW = W.algebra(K)
                sage: E = KW.demazure_lusztig_eigenvectors(q1,q2)
                sage: E.keys()
                Weyl Group of type ['B', 2] (as a matrix group acting on the ambient space)
                sage: w = W.an_element()
                sage: E[w]
                (q2/(-q1+q2))*2121 + ((-q2)/(-q1+q2))*121 - 212 + 12
            """
            W = self.basis().keys()
            if not W.cartan_type().is_finite():
                raise ValueError("the Demazure-Lusztig eigenvectors are only defined for finite Coxeter groups")
            result = self.demazure_lusztig_operators(q1, q2, affine=True).Y_eigenvectors()
            w0 = W.long_element()
            result.affine_lift = w0._mul_
            result.affine_retract = w0._mul_
            return result

