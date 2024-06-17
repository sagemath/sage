# sage.doctest: needs sage.combinat sage.modules
"""
Generic dual bases symmetric functions
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2012 Mike Zabrocki <mike.zabrocki@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.matrix.constructor import matrix
import sage.combinat.partition
import sage.data_structures.blas_dict as blas
from . import classical


class SymmetricFunctionAlgebra_dual(classical.SymmetricFunctionAlgebra_classical):
    @staticmethod
    def __classcall__(cls, dual_basis, scalar, scalar_name='', basis_name=None, prefix=None):
        """
        Normalize the arguments.

        TESTS::

            sage: w = SymmetricFunctions(QQ).w()
            sage: B1 = w.dual_basis()
            sage: B2 = w.dual_basis(prefix='d_w')
            sage: B1 is B2
            True
        """
        if prefix is None:
            prefix = 'd_'+dual_basis.prefix()
        return super().__classcall__(cls, dual_basis, scalar, scalar_name, basis_name, prefix)

    def __init__(self, dual_basis, scalar, scalar_name, basis_name, prefix):
        r"""
        Generic dual basis of a basis of symmetric functions.

        INPUT:

        - ``dual_basis`` -- a basis of the ring of symmetric functions

        - ``scalar`` -- a function `z` on partitions which determines the
          scalar product on the power sum basis by
          `\langle p_{\mu}, p_{\mu} \rangle = z(\mu)`. (Independently on the
          function chosen, the power sum basis will always be orthogonal; the
          function ``scalar`` only determines the norms of the basis elements.)
          This defaults to the function ``zee`` defined in
          ``sage.combinat.sf.sfa``, that is, the function is defined by:

          .. MATH::

              \lambda \mapsto \prod_{i = 1}^\infty m_i(\lambda)!
              i^{m_i(\lambda)}`,

          where `m_i(\lambda)` means the number of times `i` appears in
          `\lambda`. This default function gives the standard Hall scalar
          product on the ring of symmetric functions.

        - ``scalar_name`` -- (default: the empty string) a string giving a
          description of the scalar product specified by the parameter
          ``scalar``

        - ``basis_name`` -- (optional) a string to serve as name for the basis
          to be generated (such as "forgotten" in "the forgotten basis"); don't
          set it to any of the already existing basis names (such as
          ``homogeneous``, ``monomial``, ``forgotten``, etc.).

        - ``prefix`` -- (default: ``'d'`` and the prefix for ``dual_basis``)
          a string to use as the symbol for the basis

        OUTPUT:

        The basis of the ring of symmetric functions dual to the basis
        ``dual_basis`` with respect to the scalar product determined
        by ``scalar``.

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ).e()
            sage: f = e.dual_basis(prefix='m', basis_name="Forgotten symmetric functions"); f
            Symmetric Functions over Rational Field in the Forgotten symmetric functions basis
            sage: TestSuite(f).run(elements=[f[1,1]+2*f[2], f[1]+3*f[1,1]])
            sage: TestSuite(f).run() # long time (11s on sage.math, 2011)

        This class defines canonical coercions between ``self`` and
        ``self^*``, as follow:

        Lookup for the canonical isomorphism from ``self`` to `P`
        (=powersum), and build the adjoint isomorphism from `P^*` to
        ``self^*``. Since `P` is self-adjoint for this scalar product,
        derive an isomorphism from `P` to ``self^*``, and by composition
        with the above get an isomorphism from ``self`` to ``self^*`` (and
        similarly for the isomorphism ``self^*`` to ``self``).

        This should be striped down to just (auto?) defining canonical
        isomorphism by adjunction (as in MuPAD-Combinat), and let
        the coercion handle the rest.

        Inversions may not be possible if the base ring is not a field::

            sage: m = SymmetricFunctions(ZZ).m()
            sage: h = m.dual_basis(lambda x: 1)
            sage: h[2,1]
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer

        By transitivity, this defines indirect coercions to and from all other bases::

            sage: s = SymmetricFunctions(QQ['t'].fraction_field()).s()
            sage: t = QQ['t'].fraction_field().gen()
            sage: zee_hl = lambda x: x.centralizer_size(t=t)
            sage: S = s.dual_basis(zee_hl)
            sage: S(s([2,1]))
            (-t/(t^5-2*t^4+t^3-t^2+2*t-1))*d_s[1, 1, 1] + ((-t^2-1)/(t^5-2*t^4+t^3-t^2+2*t-1))*d_s[2, 1] + (-t/(t^5-2*t^4+t^3-t^2+2*t-1))*d_s[3]

        TESTS:

        Regression test for :issue:`12489`. This issue improving
        equality test revealed that the conversion back from the dual
        basis did not strip cancelled terms from the dictionary::

            sage: y = e[1, 1, 1, 1] - 2*e[2, 1, 1] + e[2, 2]
            sage: sorted(f.element_class(f, dual = y))
            [([1, 1, 1, 1], 6), ([2, 1, 1], 2), ([2, 2], 1)]
        """
        self._dual_basis = dual_basis
        self._scalar = scalar
        self._scalar_name = scalar_name

        # Set up the cache

        # cache for the coordinates of the elements
        # of ``dual_basis`` with respect to ``self``
        self._to_self_cache = {}
        # cache for the coordinates of the elements
        # of ``self`` with respect to ``dual_basis``
        self._from_self_cache = {}
        # cache for transition matrices which contain the coordinates of
        # the elements of ``dual_basis`` with respect to ``self``
        self._transition_matrices = {}
        # cache for transition matrices which contain the coordinates of
        # the elements of ``self`` with respect to ``dual_basis``
        self._inverse_transition_matrices = {}

        scalar_target = scalar(sage.combinat.partition.Partition([1])).parent()
        scalar_target = (scalar_target.one()*dual_basis.base_ring().one()).parent()

        self._sym = sage.combinat.sf.sf.SymmetricFunctions(scalar_target)
        self._p = self._sym.power()

        classical.SymmetricFunctionAlgebra_classical.__init__(self, self._sym,
                                                              basis_name=basis_name,
                                                              prefix=prefix)

        # temporary until Hom(GradedHopfAlgebrasWithBasis work better)
        category = sage.categories.all.ModulesWithBasis(self.base_ring())
        self.register_coercion(SetMorphism(Hom(self._dual_basis, self, category), self._dual_to_self))
        self._dual_basis.register_coercion(SetMorphism(Hom(self, self._dual_basis, category), self._self_to_dual))

    def construction(self):
        """
        Return a pair ``(F, R)``, where ``F`` is a
        :class:`SymmetricFunctionsFunctor` and `R` is a ring, such
        that ``F(R)`` returns ``self``.

        EXAMPLES::

            sage: w = SymmetricFunctions(ZZ).witt()
            sage: w.dual_basis().construction()
            (SymmetricFunctionsFunctor[dual Witt], Integer Ring)
        """
        return DualBasisFunctor(self), self.base_ring()

    def _dual_to_self(self, x):
        """
        Coerce an element of the dual of ``self`` canonically into ``self``.

        INPUT:

        - ``x`` -- an element in the dual basis of ``self``

        OUTPUT: ``x`` expressed in the basis ``self``

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).monomial()
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: h._dual_to_self(m([2,1]) + 3*m[1,1,1])
            d_m[1, 1, 1] - d_m[2, 1]
            sage: hh = m.realization_of().h()
            sage: h._dual_to_self(m(hh([2,2,2])))
            d_m[2, 2, 2]

        :: Note that the result is not correct if ``x`` is not an element of the
        dual basis of ``self``

            sage: h._dual_to_self(m([2,1]))
            -2*d_m[1, 1, 1] + 5*d_m[2, 1] - 3*d_m[3]
            sage: h._dual_to_self(hh([2,1]))
            -2*d_m[1, 1, 1] + 5*d_m[2, 1] - 3*d_m[3]

        This is for internal use only. Please use instead::

            sage: h(m([2,1]) + 3*m[1,1,1])
            d_m[1, 1, 1] - d_m[2, 1]
        """
        return self._element_class(self, dual=x)

    def _self_to_dual(self, x):
        """
        Coerce an element of ``self`` canonically into the dual.

        INPUT:

        - ``x`` -- an element of ``self``

        OUTPUT: ``x`` expressed in the dual basis

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).monomial()
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: h._self_to_dual(h([2,1]) + 3*h[1,1,1])
            21*m[1, 1, 1] + 11*m[2, 1] + 4*m[3]

        This is for internal use only. Please use instead::

            sage: m(h([2,1]) + 3*h[1,1,1])
            21*m[1, 1, 1] + 11*m[2, 1] + 4*m[3]

        or::

            sage: (h([2,1]) + 3*h[1,1,1]).dual()
            21*m[1, 1, 1] + 11*m[2, 1] + 4*m[3]
        """
        return x.dual()

    def _dual_basis_default(self):
        """
        Return the default value for ``self.dual_basis()``.

        This returns the basis ``self`` has been built from by
        duality.

        .. WARNING::

            This is not necessarily the dual basis for the standard
            (Hall) scalar product!

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).monomial()
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: h.dual_basis()
            Symmetric Functions over Rational Field in the monomial basis
            sage: m2 = h.dual_basis(zee, prefix='m2')
            sage: m([2])^2
            2*m[2, 2] + m[4]
            sage: m2([2])^2
            2*m2[2, 2] + m2[4]

        TESTS::

            sage: h.dual_basis() is h._dual_basis_default()
            True
        """
        return self._dual_basis

    def basis_name(self):
        r"""
        Return the name of the basis of ``self``.

        This is used for output and, for the classical bases of
        symmetric functions, to connect this basis with :ref:`Symmetrica <spkg_symmetrica>`.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: f = Sym.f()
            sage: f.basis_name()
            'forgotten'
        """
        if self._basis_name is None:
            return "dual " + self._dual_basis.basis_name()
        return self._basis_name

    def _repr_(self):
        """
        Representation of ``self``.

        OUTPUT: string description of ``self``

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).monomial()
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee); h #indirect doctests
            Dual basis to Symmetric Functions over Rational Field in the monomial basis
            sage: h = m.dual_basis(scalar=zee, scalar_name='Hall scalar product'); h #indirect doctest
            Dual basis to Symmetric Functions over Rational Field in the monomial basis with respect to the Hall scalar product
        """
        if self._basis_name is not None:
            return super()._repr_()
        if self._scalar_name:
            return "Dual basis to %s" % self._dual_basis + " with respect to the " + self._scalar_name
        return "Dual basis to %s" % self._dual_basis

    def _precompute(self, n):
        """
        Compute the transition matrices between ``self`` and its dual basis for
        the homogeneous component of size `n`. The result is not returned,
        but stored in the cache.

        INPUT:

        - ``n`` -- nonnegative integer

        EXAMPLES::

            sage: e = SymmetricFunctions(QQ['t']).elementary()
            sage: f = e.dual_basis()
            sage: f._precompute(0)
            sage: f._precompute(1)
            sage: f._precompute(2)
            sage: l = lambda c: [ (i[0],[j for j in sorted(i[1].items())]) for i in sorted(c.items())]
            sage: l(f._to_self_cache) # note: this may depend on possible previous computations!
            [([], [([], 1)]), ([1], [([1], 1)]), ([1, 1], [([1, 1], 2), ([2], 1)]), ([2], [([1, 1], 1), ([2], 1)])]
            sage: l(f._from_self_cache)
            [([], [([], 1)]), ([1], [([1], 1)]), ([1, 1], [([1, 1], 1), ([2], -1)]), ([2], [([1, 1], -1), ([2], 2)])]
            sage: f._transition_matrices[2]
            [1 1]
            [1 2]
            sage: f._inverse_transition_matrices[2]
            [ 2 -1]
            [-1  1]
        """
        base_ring = self.base_ring()
        zero = base_ring.zero()

        # Handle the n == 0 and n == 1 cases separately
        if n == 0 or n == 1:
            part = sage.combinat.partition.Partition([1]*n)
            self._to_self_cache[ part ] = { part: base_ring.one() }
            self._from_self_cache[ part ] = { part: base_ring.one() }
            self._transition_matrices[n] = matrix(base_ring, [[1]])
            self._inverse_transition_matrices[n] = matrix(base_ring, [[1]])
            return

        partitions_n = sage.combinat.partition.Partitions_n(n).list()

        # We now get separated into two cases, depending on whether we can
        # use the power-sum basis to compute the matrix, or we have to use
        # the Schur basis.

        from sage.rings.rational_field import RationalField
        if (not base_ring.has_coerce_map_from(RationalField())) and self._scalar == sage.combinat.sf.sfa.zee:
            # This is the case when (due to the base ring not being a
            # \QQ-algebra) we cannot use the power-sum basis,
            # but (due to zee being the standard zee function) we can
            # use the Schur basis.

            schur = self._sym.schur()

            # Get all the basis elements of the n^th homogeneous component
            # of the dual basis and express them in the Schur basis
            d = {}
            for part in partitions_n:
                d[part] = schur(self._dual_basis(part))._monomial_coefficients

            # This contains the data for the transition matrix from the
            # dual basis to self.
            transition_matrix_n = matrix(base_ring, len(partitions_n), len(partitions_n))

            # This first section calculates how the basis elements of the
            # dual basis are expressed in terms of self's basis.

            # For every partition p of size n, compute self(p) in
            # terms of the dual basis using the scalar product.
            i = 0
            for s_part in partitions_n:
                # s_part corresponds to self(dual_basis(part))
                # s_mcs  corresponds to self(dual_basis(part))._monomial_coefficients
                s_mcs = {}

                # We need to compute the scalar product of d[s_part] and
                # all of the d[p_part]'s
                j = 0
                for p_part in partitions_n:
                    # Compute the scalar product of d[s_part] and d[p_part]
                    sp = zero
                    for ds_part in d[s_part]:
                        if ds_part in d[p_part]:
                            sp += d[s_part][ds_part]*d[p_part][ds_part]
                    if sp != zero:
                        s_mcs[p_part] = sp
                        transition_matrix_n[i,j] = sp

                    j += 1

                self._to_self_cache[ s_part ] = s_mcs
                i += 1

        else:
            # Now the other case. Note that just being in this case doesn't
            # guarantee that we can use the power-sum basis, but we can at
            # least try.

            # Get all the basis elements of the n^th homogeneous component
            # of the dual basis and express them in the power-sum basis
            d = {}
            for part in partitions_n:
                d[part] = self._p(self._dual_basis(part))._monomial_coefficients

            # This contains the data for the transition matrix from the
            # dual basis to self.
            transition_matrix_n = matrix(base_ring, len(partitions_n), len(partitions_n))

            # This first section calculates how the basis elements of the
            # dual basis are expressed in terms of self's basis.

            # For every partition p of size n, compute self(p) in
            # terms of the dual basis using the scalar product.
            i = 0
            for s_part in partitions_n:
                # s_part corresponds to self(dual_basis(part))
                # s_mcs  corresponds to self(dual_basis(part))._monomial_coefficients
                s_mcs = {}

                # We need to compute the scalar product of d[s_part] and
                # all of the d[p_part]'s
                j = 0
                for p_part in partitions_n:
                    # Compute the scalar product of d[s_part] and d[p_part]
                    sp = zero
                    for ds_part in d[s_part]:
                        if ds_part in d[p_part]:
                            sp += d[s_part][ds_part]*d[p_part][ds_part]*self._scalar(ds_part)
                    if sp != zero:
                        s_mcs[p_part] = sp
                        transition_matrix_n[i,j] = sp

                    j += 1

                self._to_self_cache[ s_part ] = s_mcs
                i += 1

        # Save the transition matrix
        self._transition_matrices[n] = transition_matrix_n

        # This second section calculates how the basis elements of
        # self expand in terms of the dual basis.  We do this by
        # computing the inverse of the matrix obtained above.
        inverse_transition = ~transition_matrix_n

        for i in range(len(partitions_n)):
            d_mcs = {}
            for j in range(len(partitions_n)):
                if inverse_transition[i,j] != zero:
                    d_mcs[ partitions_n[j] ] = inverse_transition[i,j]

            self._from_self_cache[ partitions_n[i] ] = d_mcs

        self._inverse_transition_matrices[n] = inverse_transition

    def transition_matrix(self, basis, n):
        r"""
        Return the transition matrix between the `n`-th homogeneous components
        of ``self`` and ``basis``.

        INPUT:

        - ``basis`` -- a target basis of the ring of symmetric functions
        - ``n`` -- nonnegative integer

        OUTPUT:

        - A transition matrix from ``self`` to ``basis`` for the elements
          of degree ``n``.  The indexing order of the rows and
          columns is the order of ``Partitions(n)``.

        EXAMPLES::

            sage: Sym = SymmetricFunctions(QQ)
            sage: s = Sym.schur()
            sage: e = Sym.elementary()
            sage: f = e.dual_basis()
            sage: f.transition_matrix(s, 5)
            [ 1 -1  0  1  0 -1  1]
            [-2  1  1 -1 -1  1  0]
            [-2  2 -1 -1  1  0  0]
            [ 3 -1 -1  1  0  0  0]
            [ 3 -2  1  0  0  0  0]
            [-4  1  0  0  0  0  0]
            [ 1  0  0  0  0  0  0]
            sage: Partitions(5).list()
            [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
            sage: s(f[2,2,1])
            s[3, 2] - 2*s[4, 1] + 3*s[5]
            sage: e.transition_matrix(s, 5).inverse().transpose()
            [ 1 -1  0  1  0 -1  1]
            [-2  1  1 -1 -1  1  0]
            [-2  2 -1 -1  1  0  0]
            [ 3 -1 -1  1  0  0  0]
            [ 3 -2  1  0  0  0  0]
            [-4  1  0  0  0  0  0]
            [ 1  0  0  0  0  0  0]
        """
        if n not in self._transition_matrices:
            self._precompute(n)

        if basis is self._dual_basis:
            return self._inverse_transition_matrices[n]
        else:
            return self._inverse_transition_matrices[n]*self._dual_basis.transition_matrix(basis, n)

    def product(self, left, right):
        """
        Return product of ``left`` and ``right``.

        Multiplication is done by performing the multiplication in the dual
        basis of ``self`` and then converting back to ``self``.

        INPUT:

        - ``left``, ``right`` -- elements of ``self``

        OUTPUT: the product of ``left`` and ``right`` in the basis ``self``

        EXAMPLES::

            sage: m = SymmetricFunctions(QQ).monomial()
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: a = h([2])
            sage: b = a*a; b # indirect doctest
            d_m[2, 2]
            sage: b.dual()
            6*m[1, 1, 1, 1] + 4*m[2, 1, 1] + 3*m[2, 2] + 2*m[3, 1] + m[4]
        """
        # Do the multiplication in the dual basis
        # and then convert back to self.
        eclass = left.__class__
        d_product = left.dual() * right.dual()
        return eclass(self, dual=d_product)

    class Element(classical.SymmetricFunctionAlgebra_classical.Element):
        """
        An element in the dual basis.

        INPUT:

        At least one of the following must be specified. The one (if
        any) which is not provided will be computed.

        - ``dictionary`` -- an internal dictionary for the
          monomials and coefficients of ``self``

        - ``dual`` -- self as an element of the dual basis
        """

        def __init__(self, A, dictionary=None, dual=None):
            """
            Create an element of a dual basis.

            TESTS::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(scalar=zee, prefix='h')
                sage: a = h([2])
                sage: ec = h._element_class
                sage: ec(h, dual=m([2]))
                -h[1, 1] + 2*h[2]
                sage: h(m([2]))
                -h[1, 1] + 2*h[2]
                sage: h([2])
                h[2]
                sage: h([2])._dual
                m[1, 1] + m[2]
                sage: m(h([2]))
                m[1, 1] + m[2]
            """
            if dictionary is None and dual is None:
                raise ValueError("you must specify either x or dual")

            parent = A
            base_ring = parent.base_ring()
            zero = base_ring.zero()

            if dual is None:
                # We need to compute the dual
                dual_dict = {}
                from_self_cache = parent._from_self_cache

                # Get the underlying dictionary for self
                s_mcs = dictionary

                # Make sure all the conversions from self to
                # to the dual basis have been precomputed
                for part in s_mcs:
                    if part not in from_self_cache:
                        parent._precompute(sum(part))

                # Create the monomial coefficient dictionary from the
                # the monomial coefficient dictionary of dual
                for s_part in s_mcs:
                    from_dictionary = from_self_cache[s_part]
                    for part in from_dictionary:
                        dual_dict[ part ] = dual_dict.get(part, zero) + base_ring(s_mcs[s_part]*from_dictionary[part])

                dual = parent._dual_basis._from_dict(dual_dict)

            if dictionary is None:
                # We need to compute the monomial coefficients dictionary
                dictionary = {}
                to_self_cache = parent._to_self_cache

                # Get the underlying dictionary for the dual
                d_mcs = dual._monomial_coefficients

                # Make sure all the conversions from the dual basis
                # to self have been precomputed
                for part in d_mcs:
                    if part not in to_self_cache:
                        parent._precompute(sum(part))

                # Create the monomial coefficient dictionary from the
                # the monomial coefficient dictionary of dual
                dictionary = blas.linear_combination( (to_self_cache[d_part], d_mcs[d_part]) for d_part in d_mcs)

            # Initialize self
            self._dual = dual
            classical.SymmetricFunctionAlgebra_classical.Element.__init__(self, A, dictionary)

        def dual(self):
            """
            Return ``self`` in the dual basis.

            OUTPUT:

            - the element ``self`` expanded in the dual basis to ``self.parent()``

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(scalar=zee)
                sage: a = h([2,1])
                sage: a.parent()
                Dual basis to Symmetric Functions over Rational Field in the monomial basis
                sage: a.dual()
                3*m[1, 1, 1] + 2*m[2, 1] + m[3]
            """
            return self._dual

        def omega(self):
            r"""
            Return the image of ``self`` under the omega automorphism.

            The *omega automorphism* is defined to be the unique algebra
            endomorphism `\omega` of the ring of symmetric functions that
            satisfies `\omega(e_k) = h_k` for all positive integers `k`
            (where `e_k` stands for the `k`-th elementary symmetric
            function, and `h_k` stands for the `k`-th complete homogeneous
            symmetric function). It furthermore is a Hopf algebra
            endomorphism and an involution, and it is also known as the
            *omega involution*. It sends the power-sum symmetric function
            `p_k` to `(-1)^{k-1} p_k` for every positive integer `k`.

            The images of some bases under the omega automorphism are given by

            .. MATH::

                \omega(e_{\lambda}) = h_{\lambda}, \qquad
                \omega(h_{\lambda}) = e_{\lambda}, \qquad
                \omega(p_{\lambda}) = (-1)^{|\lambda| - \ell(\lambda)}
                p_{\lambda}, \qquad
                \omega(s_{\lambda}) = s_{\lambda^{\prime}},

            where `\lambda` is any partition, where `\ell(\lambda)` denotes
            the length (:meth:`~sage.combinat.partition.Partition.length`)
            of the partition `\lambda`, where `\lambda^{\prime}` denotes the
            conjugate partition
            (:meth:`~sage.combinat.partition.Partition.conjugate`) of
            `\lambda`, and where the usual notations for bases are used
            (`e` = elementary, `h` = complete homogeneous, `p` = powersum,
            `s` = Schur).

            :meth:`omega_involution` is a synonym for the :meth:`omega`
            method.

            OUTPUT: the result of applying omega to ``self``

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(zee)
                sage: hh = SymmetricFunctions(QQ).homogeneous()
                sage: hh([2,1]).omega()
                h[1, 1, 1] - h[2, 1]
                sage: h([2,1]).omega()
                d_m[1, 1, 1] - d_m[2, 1]
            """
            eclass = self.__class__
            return eclass(self.parent(), dual=self._dual.omega() )

        omega_involution = omega

        def scalar(self, x):
            """
            Return the standard scalar product of ``self`` and ``x``.

            INPUT:

            - ``x`` -- element of the symmetric functions

            OUTPUT: the scalar product between ``x`` and ``self``

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(scalar=zee)
                sage: a = h([2,1])
                sage: a.scalar(a)
                2
            """
            return self._dual.scalar(x)

        def scalar_hl(self, x):
            """
            Return the Hall-Littlewood scalar product of ``self`` and ``x``.

            INPUT:

            - ``x`` -- element of the same dual basis as ``self``

            OUTPUT: the Hall-Littlewood scalar product between ``x`` and ``self``

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(scalar=zee)
                sage: a = h([2,1])
                sage: a.scalar_hl(a)
                (-t - 2)/(t^4 - 2*t^3 + 2*t - 1)
            """
            return self._dual.scalar_hl(x)

        def _add_(self, y):
            """
            Add two elements in the dual basis.

            INPUT:

            - ``y`` -- element of the same dual basis as ``self``

            OUTPUT: the sum of ``self`` and ``y``

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(zee)
                sage: a = h([2,1])+h([3]); a # indirect doctest
                d_m[2, 1] + d_m[3]
                sage: h[2,1]._add_(h[3])
                d_m[2, 1] + d_m[3]
                sage: a.dual()
                4*m[1, 1, 1] + 3*m[2, 1] + 2*m[3]
            """
            eclass = self.__class__
            return eclass(self.parent(), dual=(self.dual()+y.dual()))

        def _neg_(self):
            """
            Return the negative of ``self``.

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(zee)
                sage: -h([2,1]) # indirect doctest
                -d_m[2, 1]
            """
            eclass = self.__class__
            return eclass(self.parent(), dual=self.dual()._neg_())

        def _sub_(self, y):
            """
            Subtract two elements in the dual basis.

            INPUT:

            - ``y`` -- element of the same dual basis as ``self``

            OUTPUT: the difference of ``self`` and ``y``

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(zee)
                sage: h([2,1])-h([3]) # indirect doctest
                d_m[2, 1] - d_m[3]
                sage: h[2,1]._sub_(h[3])
                d_m[2, 1] - d_m[3]
            """
            eclass = self.__class__
            return eclass(self.parent(), dual=(self.dual()-y.dual()))

        def _div_(self, y):
            """
            Divide an element ``self`` of the dual basis by ``y``.

            INPUT:

            - ``y`` -- element of base field

            OUTPUT: the element ``self`` divided by ``y``

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(zee)
                sage: a = h([2,1])+h([3])
                sage: a/2 # indirect doctest
                1/2*d_m[2, 1] + 1/2*d_m[3]
            """
            return self*(~y)

        def __invert__(self):
            """
            Invert ``self`` (only possible if ``self`` is a scalar
            multiple of `1` and we are working over a field).

            OUTPUT: multiplicative inverse of ``self`` if possible

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(zee)
                sage: a = h(2); a
                2*d_m[]
                sage: ~a
                1/2*d_m[]
                sage: a = 3*h[1]
                sage: a.__invert__()
                Traceback (most recent call last):
                ...
                ValueError: cannot invert self (= 3*m[1])
            """
            eclass = self.__class__
            return eclass(self.parent(), dual=~self.dual())

        def expand(self, n, alphabet='x'):
            """
            Expand the symmetric function ``self`` as a symmetric polynomial
            in ``n`` variables.

            INPUT:

            - ``n`` -- nonnegative integer

            - ``alphabet`` -- (default: ``'x'``) a variable for the expansion

            OUTPUT:

            A monomial expansion of ``self`` in the `n` variables
            labelled by ``alphabet``.

            EXAMPLES::

                sage: m = SymmetricFunctions(QQ).monomial()
                sage: zee = sage.combinat.sf.sfa.zee
                sage: h = m.dual_basis(zee)
                sage: a = h([2,1])+h([3])
                sage: a.expand(2)
                2*x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + 2*x1^3
                sage: a.dual().expand(2)
                2*x0^3 + 3*x0^2*x1 + 3*x0*x1^2 + 2*x1^3
                sage: a.expand(2,alphabet='y')
                2*y0^3 + 3*y0^2*y1 + 3*y0*y1^2 + 2*y1^3
                sage: a.expand(2,alphabet='x,y')
                2*x^3 + 3*x^2*y + 3*x*y^2 + 2*y^3
                sage: h([1]).expand(0)
                0
                sage: (3*h([])).expand(0)
                3
            """
            return self._dual.expand(n, alphabet)


from sage.combinat.sf.sfa import SymmetricFunctionsFunctor
class DualBasisFunctor(SymmetricFunctionsFunctor):
    """
    A constructor for algebras of symmetric functions constructed by
    duality.

    EXAMPLES::

        sage: w = SymmetricFunctions(ZZ).witt()
        sage: w.dual_basis().construction()
        (SymmetricFunctionsFunctor[dual Witt], Integer Ring)
    """
    def __init__(self, basis):
        r"""
        Initialize the functor.

        INPUT:

        - ``basis`` -- the basis of the symmetric function algebra

        TESTS::

            sage: w = SymmetricFunctions(ZZ).witt()
            sage: F = w.dual_basis().construction()[0]
            sage: TestSuite(F).run()
        """
        self._dual_basis = basis._dual_basis
        self._basis_name = basis._basis_name
        self._scalar = basis._scalar
        self._scalar_name = basis._scalar_name
        self._prefix = basis._prefix
        super().__init__(basis, self._basis_name)

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        EXAMPLES::

            sage: m = SymmetricFunctions(ZZ).monomial()
            sage: zee = sage.combinat.sf.sfa.zee
            sage: h = m.dual_basis(scalar=zee)
            sage: F, R = h.construction()  # indirect doctest
            sage: F(QQ)
            Dual basis to Symmetric Functions over Rational Field in the monomial basis

            sage: b = m.dual_basis(scalar=zee).dual_basis(scalar=lambda x: 1)
            sage: F, R = b.construction()  # indirect doctest
            sage: F(QQ)
            Dual basis to Dual basis to Symmetric Functions over Rational Field in the monomial basis
        """
        dual_basis = self._dual_basis.change_ring(R)
        return self._basis(dual_basis, self._scalar, self._scalar_name,
                           self._basis_name, self._prefix)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: w = SymmetricFunctions(ZZ).witt()
            sage: w.dual_basis().construction()
            (SymmetricFunctionsFunctor[dual Witt], Integer Ring)
        """
        if self._basis_name is None:
            name = "dual " + self._dual_basis.basis_name()
        else:
            name = self._basis_name
        return "SymmetricFunctionsFunctor[" + name + "]"


# Backward compatibility for unpickling
from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.sf.dual', 'SymmetricFunctionAlgebraElement_dual',  SymmetricFunctionAlgebra_dual.Element)
