"""
Fusion Rings
"""
# ****************************************************************************
#  Copyright (C) 2023 Guillermo Aboumrad <gh_willieab>
#                     Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import product, zip_longest
from multiprocessing import Pool, set_start_method
from sage.combinat.q_analogues import q_int
from sage.algebras.fusion_rings.fast_parallel_fusion_ring_braid_repn import (
    executor,
    _unflatten_entries
)
from sage.categories.all import Algebras, AlgebrasWithBasis, Groups
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.weyl_characters import WeylCharacterRing
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.misc import inject_variable
from sage.rings.integer_ring import ZZ
from sage.rings.number_field.number_field import CyclotomicField
from sage.rings.qqbar import QQbar

class FusionRing(CombinatorialFreeModule):
    @staticmethod
    def __classcall_private__(cls, input_data, level=None, base_ring=ZZ, prefix=None, conjugate=False, cyclotomic_order=None, fusion_labels=None, inject_variables=False, **kwds):
        """
        Select the correct parent depending on the given input.
        """
        if input_data in Groups:
            if prefix is None:
                prefix = "s"
            from sage.algebras.fusion_rings.fusion_double import FusionDouble
            return FusionDouble(input_data, base_ring=base_ring, prefix=prefix, cyclotomic_order=cyclotomic_order, fusion_labels=fusion_labels, inject_variables=inject_variables)
        else:
            try:
                ct = CartanType(input_data)
            except ValueError:
                raise ValueError("Input data must be a group or an object that can be coerced into a CartanType")
            assert level is not None, "Fusion level must be supplied when input data is a CartanType"
            from sage.algebras.fusion_rings.verlinde_algebra import VerlindeAlgebra
            return VerlindeAlgebra(ct, level, conjugate=conjugate, base_ring=base_ring, prefix=prefix, cyclotomic_order=cyclotomic_order, fusion_labels=fusion_labels, inject_variables=inject_variables)

    def __init__(self, names, base_ring=ZZ, prefix=None, conjugate=False, cyclotomic_order=None, fusion_labels=None, inject_variables=False):
        self._names = names
        self._rank = len(names)
        self._conj = -1 if conjugate else 1
        self._cyclotomic_order = cyclotomic_order
        self._field = None
        self._basecoer = None
        cat = AlgebrasWithBasis(base_ring).Subobjects()
        CombinatorialFreeModule.__init__(self, base_ring, [k for k in self._names], category=cat)
        self._fusion_labels = fusion_labels
        self.fusion_labels(labels=fusion_labels, inject_variables=inject_variables)

    @abstract_method
    def _repr_(self):
        pass

    @abstract_method
    def Nk_ij(self, i, j, k):
        r"""
        Returns the fusion coefficient `N^k_{ij}`
        """
        pass

    @abstract_method
    def one(self):
        """
        Get the multiplicative identity of ``self``.
        """
        pass

    @abstract_method
    def s_ij(self, elt_i, elt_j, unitary=False, base_coercion=True):
        """
        Get the `(i, j)`-entry of the S-matrix of ``self``.
        """
        pass

    def _element_constructor(self, k):
        """
        Construct a monomial (basis element) from a key.

        INPUT:

        - ``key`` -- a key for the dictionary `self._names`

        EXAMPLES::

            sage: F = FusionRing(SymmetricGroup(3), prefix="n")
            sage: F._names
            {0: 'n0', 1: 'n1', 2: 'n2', 3: 'n3', 4: 'n4', 5: 'n5', 6: 'n6', 7: 'n7'}
            sage: [F._element_constructor(x) for x in F._names]
            [n0, n1, n2, n3, n4, n5, n6, n7]
        """
        return self.monomial(k)

    def _repr_term(self, t):
        """
        EXAMPLES::

            sage: F = FusionRing(CyclicPermutationGroup(2))
            sage: [F._repr_term(t) for t in F._names]
            ['s0', 's1', 's2', 's3']
        """
        if self._fusion_labels is not None:
            idx = self.get_order().index(t)
            return self._fusion_labels[idx]
        return self._names[t]

    def conj_matrix(self):
        r"""
        Return the conjugation matrix, which is the permutation matrix
        for the conjugation (dual) operation on basis elements.

        EXAMPLES::

            sage: FusionRing("A2", 1).conj_matrix()
            [1 0 0]
            [0 0 1]
            [0 1 0]
        """
        b = self.basis().list()
        return matrix(ZZ, [[i == j.dual() for i in b] for j in b])

    def D_plus(self, base_coercion=True):
        r"""
        Return `\sum d_i^2\theta_i` where `i` runs through the simple objects,
        `d_i` is the quantum dimension and `\theta_i` is the twist.

        This is denoted `p_+` in [BaKi2001]_ Chapter 3.

        EXAMPLES::

            sage: B31 = FusionRing("B3", 1)
            sage: Dp = B31.D_plus(); Dp
            2*zeta48^13 - 2*zeta48^5
            sage: Dm = B31.D_minus(); Dm
            -2*zeta48^3
            sage: Dp*Dm == B31.global_q_dimension()
            True
            sage: c = B31.virasoro_central_charge(); c
            7/2
            sage: Dp/Dm == B31.root_of_unity(c/2)
            True

        For the Drinfeld double, it equals the order of the group.

        EXAMPLES::

            sage: FusionRing(DihedralGroup(8)).D_plus()
            16
        """
        ret = sum((x.q_dimension(base_coercion=False))**2 * x.ribbon(base_coercion=False) for x in self.basis())
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def D_minus(self, base_coercion=True):
        r"""
        Return `\sum d_i^2\theta_i^{-1}` where `i` runs through the simple
        objects, `d_i` is the quantum dimension and `\theta_i` is the twist.

        This is denoted `p_-` in [BaKi2001]_ Chapter 3.

        EXAMPLES::

            sage: E83 = FusionRing("E8", 3, conjugate=True)
            sage: [Dp, Dm] = [E83.D_plus(), E83.D_minus()]
            sage: Dp*Dm == E83.global_q_dimension()
            True
            sage: c = E83.virasoro_central_charge(); c
            -248/11
            sage: Dp*Dm == E83.global_q_dimension()
            True

        For the Drinfeld double, it equals the order of the group.

        EXAMPLES::

            sage: FusionRing(DihedralGroup(9)).D_minus()
            18
        """
        ret = sum((x.q_dimension(base_coercion=False))**2 / x.ribbon(base_coercion=False) for x in self.basis())
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    @cached_method
    def dual(self, i):
        r"""
        Return the dual object ``i^\ast`` to ``i``. The dual is also
        available as an element method of ``i``.

        EXAMPLES::

            sage: K = FusionRing(CyclicPermutationGroup(3),prefix="k")
            sage: [(x,K.dual(x)) for x in K.basis()]
            [(k0, k0),
            (k1, k2),
            (k2, k1),
            (k3, k6),
            (k4, k8),
            (k5, k7),
            (k6, k3),
            (k7, k5),
            (k8, k4)]
            sage: all(K.dual(x)==x.dual() for x in K.basis())
            True
        """
        sz = self.one()
        for j in self.basis():
            if self.Nk_ij(i,j,sz) > 0:
                return j

    @cached_method
    def field(self):
        """
        Return a cyclotomic field large enough to contain the S- and R-matrices.

        EXAMPLES::

            sage: FusionRing("A2", 2).field()
            Cyclotomic Field of order 60 and degree 16
            sage: FusionRing("B2", 2).field()
            Cyclotomic Field of order 40 and degree 16
            sage: FusionRing(SymmetricGroup(3)).field()
            Cyclotomic Field of order 24 and degree 8
        """
        return CyclotomicField(4 * self._cyclotomic_order)

    def fusion_labels(self, labels=None, inject_variables=False):
        r"""
        Set the labels of the basis.

        INPUT:

        - ``labels`` -- (default: ``None``) a list of strings or string
        - ``inject_variables`` -- (default: ``False``) if ``True``, then
          inject the variable names into the global namespace; note that
          this could override objects already defined

        If ``labels`` is a list, the length of the list must equal the
        number of basis elements. These become the names of
        the basis elements.

        If ``labels`` is a string, this is treated as a prefix and a
        list of names is generated.

        If ``labels`` is ``None``, then this resets the labels to the default.

        EXAMPLES::

            sage: A13 = FusionRing("A1", 3)
            sage: A13.fusion_labels("x")
            sage: fb = list(A13.basis()); fb
            [x0, x1, x2, x3]
            sage: Matrix([[x*y for y in A13.basis()] for x in A13.basis()])
            [     x0      x1      x2      x3]
            [     x1 x0 + x2 x1 + x3      x2]
            [     x2 x1 + x3 x0 + x2      x1]
            [     x3      x2      x1      x0]

        We give an example where the variables are injected into the
        global namespace::

            sage: A13.fusion_labels("y", inject_variables=True)
            sage: y0
            y0
            sage: y0.parent() is A13
            True

        We reset the labels to the default::

            sage: A13.fusion_labels()
            sage: fb
            [A13(0), A13(1), A13(2), A13(3)]
            sage: y0
            A13(0)
        """
        # Maintain a tuple of names sorted according to (ordered) CFM basis
        if labels is None:
            # Reset the fusion labels
            labels = tuple(str(self._names[b]) for b in self.get_order())

        B = self.basis()
        if isinstance(labels, str):
            labels = tuple(labels + str(k) for k in range(len(B)))
        elif len(labels) != len(B):
            raise ValueError('invalid data')

        for j, b in enumerate(self.get_order()):
            if inject_variables:
                inject_variable(labels[j], B[b])
        self._fusion_labels = labels

    def fvars_field(self):
        r"""
        Return a field containing the ``CyclotomicField`` computed by
        :meth:`field` as well as all the F-symbols of the associated
        ``FMatrix`` factory object.

        This method is only available if ``self`` is multiplicity-free.

        OUTPUT:

        Depending on the ``CartanType`` associated to ``self`` and whether
        a call to an F-matrix solver has been made, this method
        will return the same field as :meth:`field`, a :func:`NumberField`,
        or the :class:`QQbar<AlgebraicField>`.
        See :meth:`FMatrix.attempt_number_field_computation` for more details.

        Before running an F-matrix solver, the output of this method matches
        that of :meth:`field`. However, the output may change upon successfully
        computing F-symbols. Requesting braid generators triggers a call to
        :meth:`FMatrix.find_orthogonal_solution`, so the output of this method
        may change after such a computation.

        By default, the output of methods like :meth:`r_matrix`,
        :meth:`s_matrix`, :meth:`twists_matrix`, etc. will lie in the
        ``fvars_field``, unless the ``base_coercion`` option is set to
        ``False``.

        This method does not trigger a solver run.

        EXAMPLES::

            sage: A13 = FusionRing("A1", 3, fusion_labels="a", inject_variables=True)
            sage: A13.fvars_field()
            Cyclotomic Field of order 40 and degree 16
            sage: A13.field()
            Cyclotomic Field of order 40 and degree 16
            sage: a2**4
            2*a0 + 3*a2
            sage: comp_basis, sig = A13.get_braid_generators(a2, a2, 3, verbose=False)    # long time (<3s)
            sage: A13.fvars_field()                                                       # long time
            Number Field in a with defining polynomial y^32 - ... - 500*y^2 + 25
            sage: a2.q_dimension().parent()                                               # long time
            Number Field in a with defining polynomial y^32 - ... - 500*y^2 + 25
            sage: A13.field()
            Cyclotomic Field of order 40 and degree 16

        In some cases, the :meth:`NumberField.optimized_representation()
        <sage.rings.number_field.number_field.NumberField_absolute.optimized_representation>`
        may be used to obtain a better defining polynomial for the
        computed :func:`NumberField`.
        """
        if self.is_multiplicity_free():
            return self.get_fmatrix().field()
        raise NotImplementedError("method is only available for multiplicity free fusion rings")

    def gens_satisfy_braid_gp_rels(self, sig):
        r"""
        Return ``True`` if the matrices in the list ``sig`` satisfy
        the braid relations.

        This if `n` is the cardinality of ``sig``, this
        confirms that these matrices define a representation of
        the Artin braid group on `n+1` strands. Tests correctness of
        :meth:`get_braid_generators`.

        EXAMPLES::

            sage: F41 = FusionRing("F4", 1, fusion_labels="f", inject_variables=True)
            sage: f1*f1
            f0 + f1
            sage: comp, sig = F41.get_braid_generators(f1, f0, 4, verbose=False)
            sage: F41.gens_satisfy_braid_gp_rels(sig)
            True
        """
        n = len(sig)
        braid_rels = all(sig[i] * sig[i+1] * sig[i] == sig[i+1] * sig[i] * sig[i+1] for i in range(n-1))
        far_comm = all(sig[i] * sig[j] == sig[j] * sig[i] for i, j in product(range(n), repeat=2) if abs(i-j) > 1 and i > j)
        singular = any(s.is_singular() for s in sig)
        return braid_rels and far_comm and not singular

    def get_braid_generators(self,
                            fusing_anyon,
                            total_charge_anyon,
                            n_strands,
                            checkpoint=False,
                            save_results="",
                            warm_start="",
                            use_mp=True,
                            verbose=True):
        r"""
        Compute generators of the Artin braid group on ``n_strands`` strands.

        If `a = ` ``fusing_anyon`` and `b = ` ``total_charge_anyon``
        the generators are endomorphisms of `\text{Hom}(b, a^n)`.

        INPUT:

        - ``fusing_anyon`` -- a basis element of ``self``
        - ``total_charge_anyon`` -- a basis element of ``self``
        - ``n_strands`` -- a positive integer greater than 2
        - ``checkpoint`` -- (default: ``False``) a boolean indicating
          whether the F-matrix solver should pickle checkpoints
        - ``save_results`` -- (optional) a string indicating the name of
          a file in which to pickle computed F-symbols for later use
        - ``warm_start`` -- (optional) a string indicating the name of a
          pickled checkpoint file to "warm" start the F-matrix solver.
          The pickle may be a checkpoint generated by the solver, or
          a file containing solver results. If all F-symbols are known,
          we don't run the solver again.
        - ``use_mp`` -- (default: ``True``) a boolean indicating whether
          to use multiprocessing to speed up the computation; this is
          highly recommended. Python 3.8+ is required.
        - ``verbose`` -- (default: ``True``) boolean indicating whether
          to be verbose with the computation

        For more information on the optional parameters, see
        :meth:`FMatrix.find_orthogonal_solution`.

        Given a simple object in the fusion category, here called
        ``fusing_anyon`` allowing the universal R-matrix to act on adjacent
        pairs in the fusion of ``n_strands`` copies of ``fusing_anyon``
        produces an action of the braid group. This representation can
        be decomposed over another anyon, here called ``total_charge_anyon``.
        See [CHW2015]_.

        OUTPUT:

        The method outputs a pair of data ``(comp_basis, sig)`` where
        ``comp_basis`` is a list of basis elements of the braid group
        module, parametrized by a list of fusion ring elements describing
        a fusion tree. For example with 5 strands the fusion tree
        is as follows. See :meth:`get_computational_basis`
        for more information.

        .. IMAGE:: ../../../media/fusiontree.png
           :scale: 45
           :align: center

        ``sig`` is a list of braid group generators as matrices. In
        some cases these will be represented as sparse matrices.

        In the following example we compute a 5-dimensional braid group
        representation on 5 strands associated to the spin representation
        in the modular tensor category `SU(2)_4 \cong SO(3)_2`.

        EXAMPLES::

            sage: A14 = FusionRing("A1", 4)
            sage: A14.get_order()
            [(0, 0), (1/2, -1/2), (1, -1), (3/2, -3/2), (2, -2)]
            sage: A14.fusion_labels(["one", "two", "three", "four", "five"], inject_variables=True)
            sage: [A14(x) for x in A14.get_order()]
            [one, two, three, four, five]
            sage: two ** 5
            5*two + 4*four
            sage: comp_basis, sig = A14.get_braid_generators(two, two, 5, verbose=False) # long time
            sage: A14.gens_satisfy_braid_gp_rels(sig)                                 # long time
            True
            sage: len(comp_basis) == 5                                                # long time
            True
        """
        if n_strands < 3:
            raise ValueError("the number of strands must be an integer at least 3")
        # Construct associated FMatrix object and solve for F-symbols
        self.get_fmatrix()
        if self.fmats._chkpt_status < 7:
            self.fmats.find_orthogonal_solution(checkpoint=checkpoint,
                                                save_results=save_results,
                                                warm_start=warm_start,
                                                use_mp=use_mp,
                                                verbose=verbose)

        # Set multiprocessing parameters. Context can only be set once, so we try to set it
        try:
            set_start_method('fork')
        except RuntimeError:
            pass
        # Turn off multiprocessing when field is QQbar due to pickling issues introduced by PARI upgrade in trac ticket #30537
        pool = Pool() if use_mp and self.fvars_field() != QQbar else None

        # Set up computational basis and compute generators one at a time
        a, b = fusing_anyon, total_charge_anyon
        comp_basis = self.get_computational_basis(a, b, n_strands)
        d = len(comp_basis)
        if verbose:
            print("Computing an {}-dimensional representation of the Artin braid group on {} strands...".format(d, n_strands))

        # Compute diagonal odd-indexed generators using the 3j-symbols
        gens = {2*i+1: diagonal_matrix(self.r_matrix(a, a, c[i]) for c in comp_basis) for i in range(n_strands//2)}

        # Compute even-indexed generators using F-matrices
        for k in range(1, n_strands//2):
            entries = self._emap('sig_2k', (k, a, b, n_strands), pool)

            # Build cyclotomic field element objects from tuple of rationals repn
            _unflatten_entries(self, entries)
            gens[2*k] = matrix(dict(entries))

        # If n_strands is odd, we compute the final generator
        if n_strands % 2:
            entries = self._emap('odd_one_out', (a, b, n_strands), pool)

            # Build cyclotomic field element objects from tuple of rationals repn
            _unflatten_entries(self, entries)
            gens[n_strands-1] = matrix(dict(entries))

        return comp_basis, [gens[k] for k in sorted(gens)]

    def _emap(self, mapper, input_args, worker_pool=None):
        r"""
        Apply the given mapper to each element of the given input iterable
        and return the results (with no duplicates) in a list.

        INPUT:

        - ``mapper`` -- a string specifying the name of a function defined
          in the ``fast_parallel_fusion_ring_braid_repn`` module
        - ``input_args`` -- a tuple of arguments to be passed to mapper

        This method applies the mapper in parallel if a ``worker_pool``
        is provided.

        .. NOTE::

            If ``worker_pool`` is not provided, function maps and reduces on
            a single process. If ``worker_pool`` is provided, the function
            attempts to determine whether it should use multiprocessing
            based on the length of the input iterable. If it cannot determine
            the length of the input iterable then it uses multiprocessing
            with the default chunksize of `1` if chunksize is not
            explicitly provided.

        EXAMPLES::

            sage: FR = FusionRing("A1", 4)
            sage: FR.fusion_labels(['idd', 'one', 'two', 'three', 'four'], inject_variables=True)
            sage: fmats = FR.get_fmatrix()
            sage: fmats.find_orthogonal_solution(verbose=False)             # long time
            sage: len(FR._emap('sig_2k', (1, one, one, 5)))                 # long time
            13
            sage: FR = FusionRing("A1", 2)
            sage: FR.fusion_labels("a", inject_variables=True)
            sage: fmats = FR.get_fmatrix()
            sage: fmats.find_orthogonal_solution(verbose=False)
            sage: len(FR._emap('odd_one_out', (a1, a1, 7)))
            16
        """
        n_proc = worker_pool._processes if worker_pool is not None else 1
        input_iter = [(child_id, n_proc, input_args) for child_id in range(n_proc)]
        no_mp = worker_pool is None
        # Map phase
        input_iter = zip_longest([], input_iter, fillvalue=(mapper, id(self)))
        results = list()
        if no_mp:
            mapped = map(executor, input_iter)
        else:
            mapped = worker_pool.imap_unordered(executor, input_iter, chunksize=1)
        # Reduce phase
        for worker_results in mapped:
            results.extend(worker_results)
        return results

    def get_computational_basis(self, a, b, n_strands):
        r"""
        Return the so-called computational basis for `\text{Hom}(b, a^n)`.

        INPUT:

        - ``a`` -- a basis element
        - ``b`` -- another basis element
        - ``n_strands`` -- the number of strands for a braid group

        Let `n=` ``n_strands`` and let `k` be the greatest integer `\leq n/2`.
        The braid group acts on `\text{Hom}(b, a^n)`. This action
        is computed in :meth:`get_braid_generators`. This method
        returns the computational basis in the form of a list of
        fusion trees. Each tree is represented by an `(n-2)`-tuple

        .. MATH::

            (m_1, \ldots, m_k, l_1, \ldots, l_{k-2})

        such that each `m_j` is an irreducible constituent in `a \otimes a`
        and

        .. MATH::

            \begin{array}{l}
            b \in l_{k-2} \otimes m_{k}, \\
            l_{k-2} \in l_{k-3} \otimes m_{k-1}, \\
            \cdots, \\
            l_2 \in l_1 \otimes m_3, \\
            l_1 \in m_1 \otimes m_2,
            \end{array}

        where `z \in x \otimes y` means `N_{xy}^z \neq 0`.

        As a computational device when ``n_strands`` is odd, we pad the
        vector `(m_1, \ldots, m_k)` with an additional `m_{k+1}` equal to `a`.
        However, this `m_{k+1}` does *not* appear in the output of this method.

        The following example appears in Section 3.1 of [CW2015]_.

        EXAMPLES::

            sage: A14 = FusionRing("A1", 4)
            sage: A14.get_order()
            [(0, 0), (1/2, -1/2), (1, -1), (3/2, -3/2), (2, -2)]
            sage: A14.fusion_labels(["zero", "one", "two", "three", "four"], inject_variables=True)
            sage: [A14(x) for x in A14.get_order()]
            [zero, one, two, three, four]
            sage: A14.get_computational_basis(one, two, 4)
            [(two, two), (two, zero), (zero, two)]
        """
        def _get_trees(fr, top_row, root):
            if len(top_row) == 2:
                m1, m2 = top_row
                return [[]] if fr.Nk_ij(m1, m2, root) else []
            else:
                m1, m2 = top_row[:2]
                return [tuple([l, *b]) for l in fr.basis() for b in _get_trees(fr, [l]+top_row[2:], root) if fr.Nk_ij(m1, m2, l)]

        comp_basis = list()
        for top in product((a*a).monomials(), repeat=n_strands//2):
            # If the n_strands is odd, we must extend the top row by a fusing anyon
            top_row = list(top)+[a]*(n_strands%2)
            comp_basis.extend(tuple([*top, *levels]) for levels in _get_trees(self, top_row, b))
        return comp_basis

    def get_fmatrix(self, *args, **kwargs):
        r"""
        Construct an :class:`FMatrix` factory to solve the pentagon relations
        and organize the resulting F-symbols.

        We only need this attribute to compute braid group representations.

        EXAMPLES::

            sage: A15 = FusionRing("A1", 5)
            sage: A15.get_fmatrix()
            F-Matrix factory for The Fusion Ring of Type A1 and level 5 with Integer Ring coefficients

        EXAMPLES::

            sage: f = FusionRing(SymmetricGroup(3)).get_fmatrix(); f
            F-Matrix factory for The Fusion Ring of the Drinfeld Double of Symmetric group of order 3! as a permutation group
        """
        # Initialize fresh FMatrix object. Useful if you need to reset
        # FMatrix properties and there are various FusionRing objects (unique)
        # associated to same level and algebra.
        if not hasattr(self, 'fmats') or kwargs.get('new', False):
            kwargs.pop('new', None)
            from sage.algebras.fusion_rings.f_matrix import FMatrix
            self.fmats = FMatrix(self, *args, **kwargs)
        return self.fmats

    def get_order(self):
        r"""
        Return the weights of the basis vectors in a fixed order.

        You may change the order of the basis using :meth:`CombinatorialFreeModule.set_order`

        EXAMPLES::

            sage: A15 = FusionRing("A1", 5)
            sage: w = A15.get_order(); w
            [(0, 0), (1/2, -1/2), (1, -1), (3/2, -3/2), (2, -2), (5/2, -5/2)]
            sage: A15.set_order([w[k] for k in [0, 4, 1, 3, 5, 2]])
            sage: [A15(x) for x in A15.get_order()]
            [A15(0), A15(4), A15(1), A15(3), A15(5), A15(2)]
            sage: FusionRing(SymmetricGroup(4)).get_order()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

        .. WARNING::

            This duplicates :meth:`get_order` from
            :class:`CombinatorialFreeModule` except the result
            is *not* cached. Caching of
            :meth:`CombinatorialFreeModule.get_order` causes inconsistent
            results after calling :meth:`CombinatorialFreeModule.set_order`.

        """
        if self._order is None:
            self.set_order(self.basis().keys().list())
        return self._order

    def global_q_dimension(self, base_coercion=True):
        r"""
        Return `\sum d_i^2`, where the sum is over all simple objects
        and `d_i` is the quantum dimension.

        The global `q`-dimension is a positive real number.

        EXAMPLES::

            sage: FusionRing("E6", 1).global_q_dimension()
            3

        For the Drinfeld double, it is the square of the order of the underlying quantum group.

        EXAMPLES::

            sage: G = SymmetricGroup(4)
            sage: H = FusionRing(G)
            sage: H.global_q_dimension()
            576
            sage: sum(x.q_dimension()^2 for x in H.basis())
            576
        """
        ret = sum(x.q_dimension(base_coercion=False) ** 2 for x in self.basis())
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def is_multiplicity_free(self):
        r"""
        Return ``True`` if the fusion multiplicities
        :meth:`Nk_ij` are bounded by 1.

        The :class:`FMatrix` is available only for multiplicity free
        instances of :class:`FusionRing`.

        EXAMPLES::

            sage: [FusionRing(ct, k).is_multiplicity_free() for ct in ("A1", "A2", "B2", "C3") for k in (1, 2, 3)]
            [True, True, True, True, True, False, True, True, False, True, False, False]

        EXAMPLES::

            sage: FusionRing(SymmetricGroup(3)).is_multiplicity_free()
            True
            sage: FusionRing(SymmetricGroup(4)).is_multiplicity_free()
            False
        """
        return all(self.N_ijk(i, j, k) <= 1 for i in self.basis() for j in self.basis() for k in self.basis())

    def N_ijk(self, i, j, k):
        """
        The symmetric invariant of three simple objects,
        this returns the dimension of

        .. MATH::
           Hom(i \\otimes j\\otimes k, s0)

        where `s_0` is the unit element (assuming prefix='s').

        EXAMPLES::

            sage: G23 = FusionRing("G2", 3)
            sage: G23.fusion_labels("g")
            sage: b = G23.basis().list(); b
            [g0, g1, g2, g3, g4, g5]
            sage: [(x, y, z) for x in b for y in b for z in b if G23.N_ijk(x, y, z) > 1]
            [(g3, g3, g3), (g3, g3, g4), (g3, g4, g3), (g4, g3, g3)]
            sage: all(G23.N_ijk(x, y, z) == G23.N_ijk(y, z, x) for x in b for y in b for z in b)
            True
            sage: all(G23.N_ijk(x, y, z) == G23.N_ijk(y, x, z) for x in b for y in b for z in b)
            True

        EXAMPLES::

            sage: A = FusionRing(AlternatingGroup(4), prefix="a", inject_variables=True)
            sage: [A.N_ijk(a10, a11, x) for x in A.basis()]
            [0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        """
        return self.Nk_ij(i, j, self.dual(k))

    def product_on_basis(self, a, b):
        """
        Return the product of two basis elements corresponding to keys `a` and `b`.

        INPUT:

        - ``a`, ``b`` -- keys for the dictionary ``self._names`` representing simple objects

        EXAMPLES::

            sage: S = FusionRing(SymmetricGroup(3), prefix="s", inject_variables=True)
            sage: s3*s4
            s1 + s2 + s5 + s6 + s7
            sage: S._names
            {0: 's0', 1: 's1', 2: 's2', 3: 's3', 4: 's4', 5: 's5', 6: 's6', 7: 's7'}
            sage: S.product_on_basis(3,4)
            s1 + s2 + s5 + s6 + s7
        """
        d = {k.support_of_term() : self.Nk_ij(self.monomial(a),self.monomial(b),k) for k in self.basis()}
        return self._from_dict(d)

    @cached_method
    def r_matrix(self, i, j, k, base_coercion=True):
        r"""
        Return the R-matrix entry corresponding to the subobject ``k``
        in the tensor product of ``i`` with ``j``.

        .. WARNING::

            This method only gives complete information when `N_{ij}^k = 1`
            (an important special case). Tables of MTC including R-matrices
            may be found in Section 5.3 of [RoStWa2009]_ and in [Bond2007]_.

        The R-matrix is a homomorphism `i \otimes j \rightarrow j \otimes i`.
        This may be hard to describe since the object `i \otimes j`
        may be reducible. However if `k` is a simple subobject of
        `i \otimes j` it is also a subobject of `j \otimes i`. If we fix
        embeddings `k \rightarrow i \otimes j`, `k \rightarrow j \otimes i`
        we may ask for the scalar automorphism of `k` induced by the
        R-matrix. This method computes that scalar. It is possible to
        adjust the set of embeddings `k \rightarrow i \otimes j` (called
        a *gauge*) so that this scalar equals

        .. MATH::

            \pm \sqrt{\frac{ \theta_k }{ \theta_i \theta_j }}.

        If `i \neq j`, the gauge may be used to control the sign of
        the square root. But if `i = j` then we must be careful
        about the sign. These cases are computed by a formula
        of [BDGRTW2019]_, Proposition 2.3.

        EXAMPLES::

            sage: I = FusionRing("E8", 2, conjugate=True)  # Ising MTC
            sage: I.fusion_labels(["i0", "p", "s"], inject_variables=True)
            sage: I.r_matrix(s, s, i0) == I.root_of_unity(-1/8)
            True
            sage: I.r_matrix(p, p, i0)
            -1
            sage: I.r_matrix(p, s, s) == I.root_of_unity(-1/2)
            True
            sage: I.r_matrix(s, p, s) == I.root_of_unity(-1/2)
            True
            sage: I.r_matrix(s, s, p) == I.root_of_unity(3/8)
            True

        EXAMPLES::

            sage: C = FusionRing(SymmetricGroup(3), prefix="c", inject_variables=True)
            sage: c4*c5
            c3 + c4
            sage: [C.r_matrix(c4,c5,k) for k in [c3,c4]]
            [-zeta24^6, 1]
            sage: c6^2
            c0 + c1 + c6
            sage: [C.r_matrix(c6,c6,k) for k in [c0,c1,c6]]
            [zeta3, -zeta3, -zeta3 - 1]
        """
        if self.Nk_ij(i, j, k) == 0:
            return self.field().zero() if (not base_coercion) or (self._basecoer is None) else self.fvars_field().zero()
        if i != j:
            ret = self.root_of_unity((k.twist() - i.twist() - j.twist()) / 2)
        else:
            i0 = self.one()
            B = self.basis()
            ret = sum(y.ribbon(base_coercion=False)**2 / (i.ribbon(base_coercion=False) * x.ribbon(base_coercion=False)**2)
                   * self.s_ij(i0, y, base_coercion=False) * self.s_ij(i, z, base_coercion=False) * self.s_ijconj(x, z, base_coercion=False)
                   * self.s_ijconj(k, x, base_coercion=False) * self.s_ijconj(y, z, base_coercion=False) / self.s_ij(i0, z, base_coercion=False)
                   for x in B for y in B for z in B) / (self.total_q_order(base_coercion=False)**4)
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def root_of_unity(self, r, base_coercion=True):
        r"""
        Return `e^{i\pi r}` as an element of ``self.field()`` if possible.

        INPUT:

        - ``r`` -- a rational number

        EXAMPLES::

            sage: A11 = FusionRing("A1", 1)
            sage: A11.field()
            Cyclotomic Field of order 24 and degree 8
            sage: for n in [1..7]:
            ....:     try:
            ....:         print(n, A11.root_of_unity(2/n))
            ....:     except ValueError as err:
            ....:         print(n, err)
            1 1
            2 -1
            3 zeta24^4 - 1
            4 zeta24^6
            5 not a root of unity in the field
            6 zeta24^4
            7 not a root of unity in the field

        EXAMPLES::

            sage: H = FusionRing(DihedralGroup(6))
            sage: H.field()
            Cyclotomic Field of order 24 and degree 8
            sage: for n in [1..7]:
            ....:     try:
            ....:         print (n,H.root_of_unity(2/n))
            ....:     except ValueError as err:
            ....:         print (n,err)
            ....:
            1 1
            2 -1
            3 zeta24^4 - 1
            4 zeta24^6
            5 not a root of unity in the field
            6 zeta24^4
            7 not a root of unity in the field
        """
        n = 2 * r * self._cyclotomic_order
        if n not in ZZ:
            raise ValueError("not a root of unity in the field")
        ret = self.field().gen() ** n
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def s_ijconj(self, elt_i, elt_j, unitary=False, base_coercion=True):
        """
        Return the conjugate of the element of the S-matrix given by
        ``self.s_ij(elt_i, elt_j, base_coercion=base_coercion)``.

        See :meth:`s_ij`.

        EXAMPLES::

            sage: G21 = FusionRing("G2", 1)
            sage: b = G21.basis()
            sage: [G21.s_ijconj(x, y) for x in b for y in b]
            [1, -zeta60^14 + zeta60^6 + zeta60^4, -zeta60^14 + zeta60^6 + zeta60^4, -1]

        EXAMPLES::

            sage: P = FusionRing(CyclicPermutationGroup(3), prefix="p", inject_variables=True)
            sage: P.s_ij(p1,p3)
            zeta3
            sage: P.s_ijconj(p1,p3)
            -zeta3 - 1

        TESTS::

            sage: E62 = FusionRing("E6", 2)
            sage: E62.fusion_labels("e", inject_variables=True)
            sage: E62.s_ij(e8, e1).conjugate() == E62.s_ijconj(e8, e1)
            True
            sage: F41 = FusionRing("F4", 1)
            sage: fmats = F41.get_fmatrix()
            sage: fmats.find_orthogonal_solution(verbose=False)
            sage: b = F41.basis()
            sage: all(F41.s_ijconj(x, y) == F41._basecoer(F41.s_ij(x, y, base_coercion=False).conjugate()) for x in b for y in b)
            True
            sage: G22 = FusionRing("G2", 2)
            sage: fmats = G22.get_fmatrix()
            sage: fmats.find_orthogonal_solution(verbose=False)         # long time (~11 s)
            sage: b = G22.basis()                                       # long time
            sage: all(G22.s_ijconj(x, y) == fmats.field()(G22.s_ij(x, y, base_coercion=False).conjugate()) for x in b for y in b)   # long time
            True
        """
        ret = self.s_ij(elt_i, elt_j, unitary=unitary, base_coercion=False).conjugate()
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def s_matrix(self, unitary=False, base_coercion=True):
        r"""
        Return the S-matrix of this fusion ring.

        OPTIONAL:

        - ``unitary`` -- (default: ``False``) set to ``True`` to obtain
          the unitary S-matrix

        Without the ``unitary`` parameter, this is the matrix denoted
        `\widetilde{s}` in [BaKi2001]_.

        EXAMPLES::

            sage: D91 = FusionRing("D9", 1)
            sage: D91.s_matrix()
            [          1           1           1           1]
            [          1           1          -1          -1]
            [          1          -1 -zeta136^34  zeta136^34]
            [          1          -1  zeta136^34 -zeta136^34]
            sage: S = D91.s_matrix(unitary=True); S
            [            1/2             1/2             1/2             1/2]
            [            1/2             1/2            -1/2            -1/2]
            [            1/2            -1/2 -1/2*zeta136^34  1/2*zeta136^34]
            [            1/2            -1/2  1/2*zeta136^34 -1/2*zeta136^34]
            sage: S*S.conjugate()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]

        EXAMPLES::

            sage: FusionRing(SymmetricGroup(3)).s_matrix()
            [ 1  1  2  3  3  2  2  2]
            [ 1  1  2 -3 -3  2  2  2]
            [ 2  2  4  0  0 -2 -2 -2]
            [ 3 -3  0  3 -3  0  0  0]
            [ 3 -3  0 -3  3  0  0  0]
            [ 2  2 -2  0  0  4 -2 -2]
            [ 2  2 -2  0  0 -2 -2  4]
            [ 2  2 -2  0  0 -2  4 -2]
            sage: FusionRing(SymmetricGroup(3)).s_matrix(unitary=True)
            [ 1/6  1/6  1/3  1/2  1/2  1/3  1/3  1/3]
            [ 1/6  1/6  1/3 -1/2 -1/2  1/3  1/3  1/3]
            [ 1/3  1/3  2/3    0    0 -1/3 -1/3 -1/3]
            [ 1/2 -1/2    0  1/2 -1/2    0    0    0]
            [ 1/2 -1/2    0 -1/2  1/2    0    0    0]
            [ 1/3  1/3 -1/3    0    0  2/3 -1/3 -1/3]
            [ 1/3  1/3 -1/3    0    0 -1/3 -1/3  2/3]
            [ 1/3  1/3 -1/3    0    0 -1/3  2/3 -1/3]
        """
        b = self.basis()
        S = matrix([[self.s_ij(b[x], b[y], unitary=unitary, base_coercion=base_coercion)
                     for x in self.get_order()] for y in self.get_order()])
        return S

    def total_q_order(self, base_coercion=True):
        r"""
        Return the positive square root of :meth:`self.global_q_dimension()
        <global_q_dimension>` as an element of :meth:`self.field() <field>`.

        This is implemented as `D_{+}e^{-i\pi c/4}`, where `D_+` is
        :meth:`D_plus()` and `c` is :meth:`virasoro_central_charge()`.

        EXAMPLES::

            sage: F = FusionRing("G2", 1)
            sage: tqo=F.total_q_order(); tqo
            zeta60^15 - zeta60^11 - zeta60^9 + 2*zeta60^3 + zeta60
            sage: tqo.is_real_positive()
            True
            sage: tqo^2 == F.global_q_dimension()
            True

        For the Drinfeld double of a finite group `G`, this equals the
        cardinality of `G`.

        EXAMPLES::

            sage: FusionRing(DihedralGroup(7)).total_q_order()
            14
        """
        roots = self.global_q_dimension(base_coercion=False).sqrt(extend=False, all=True)
        ret = [root for root in roots if root.is_real_positive()][0]
        if (not base_coercion) or (self._basecoer is None):
            return ret
        return self._basecoer(ret)

    def test_braid_representation(self, max_strands=6, anyon=None):
        """
        Check that we can compute valid braid group representations.

        INPUT:

        - ``max_strands`` -- (default: 6): maximum number of braid group strands
        - ``anyon`` -- (optional) run this test on this particular simple object

        Create a braid group representation using :meth:`get_braid_generators`
        and confirms the braid relations.  This test indirectly partially
        verifies the correctness of the orthogonal F-matrix solver. If the
        code were incorrect the method would not be deterministic because the
        fusing anyon is chosen randomly. (A different choice is made for each
        number of strands tested.) However the doctest is deterministic since
        it will always return ``True``. If the anyon parameter is omitted,
        a random anyon is tested for each number of strands up to ``max_strands``.

        EXAMPLES::

            sage: A21 = FusionRing("A2", 1)
            sage: A21.test_braid_representation(max_strands=4)
            True
            sage: F41 = FusionRing("F4", 1)            # long time
            sage: F41.test_braid_representation()        # long time
            True
        """
        if not self.is_multiplicity_free(): # Braid group representation is not available if self is not multiplicity free
            raise NotImplementedError("only implemented for multiplicity free fusion rings")
        b = self.basis()
        results = []
        # Test with different numbers of strands
        for n_strands in range(3, max_strands+1):
            # Randomly select a fusing anyon. Skip the identity element, since
            # its braiding matrices are trivial
            if anyon is not None:
                a = anyon
            else:
                while True:
                    a = b.random_element()
                    if a != self.one():
                        break
            pow = a ** n_strands
            d = pow.monomials()[0]
            # Try to find 'interesting' braid group reps i.e. skip 1-d reps
            for k, v in pow.monomial_coefficients().items():
                if v > 1:
                    d = self(k)
                break
            comp_basis, sig = self.get_braid_generators(a, d, n_strands, verbose=False)
            results.append(len(comp_basis) > 0)
            results.append(self.gens_satisfy_braid_gp_rels(sig))
        return all(results)

    def twists_matrix(self):
        r"""
        Return a diagonal matrix describing the twist corresponding to
        each simple object in the ``FusionRing``.

        EXAMPLES::

            sage: B21 = FusionRing("B2", 1)
            sage: [x.twist() for x in B21.basis().list()]
            [0, 1, 5/8]
            sage: [B21.root_of_unity(x.twist()) for x in B21.basis().list()]
            [1, -1, zeta32^10]
            sage: B21.twists_matrix()
            [        1         0         0]
            [        0        -1         0]
            [        0         0 zeta32^10]
        """
        B = self.basis()
        return diagonal_matrix(B[x].ribbon() for x in self.get_order())

    class Element(CombinatorialFreeModule.Element):

        @abstract_method
        def ribbon(self, base_coercion=True):
            """
            The twist or ribbon of the simple object.
            """
            pass

        def dual(self):
            """
            The dual element under the conjugation involution.

            EXAMPLES::

                sage: G = CyclicPermutationGroup(4)
                sage: H = FusionRing(G, prefix="j")
                sage: [x for x in H.basis() if x==x.dual()]
                [j0, j1, j8, j9]
            """
            if not self.is_simple_object():
                raise ValueError("Dual is only available for simple objects of a FusionRing")
            return self.parent().dual(self)

        def is_simple_object(self):
            """
            Determine whether ``self`` is a simple object (basis element) of the fusion ring.

            EXAMPLES::

                sage: A22 = FusionRing("A2", 2)
                sage: x = A22(1, 0); x
                A22(1,0)
                sage: x.is_simple_object()
                True
                sage: x^2
                A22(0,1) + A22(2,0)
                sage: (x^2).is_simple_object()
                False

            EXAMPLES::

                sage: H = FusionRing(CyclicPermutationGroup(2), prefix="g", inject_variables=True)
                sage: [x.is_simple_object() for x in [g0, g1, g0+g1]]
                [True, True, False]
            """
            return self in self.parent().basis()

        @cached_method
        def q_dimension(self, base_coercion=True):
            """
            Compute the quantum dimension as an element of the parent's base
            cyclotomic field.

            EXAMPLES::

                sage: B22 = FusionRing("B2", 2)
                sage: [(b.q_dimension())^2 for b in B22.basis()]
                [1, 4, 5, 1, 5, 4]

            EXAMPLES::

                sage: G = AlternatingGroup(4)
                sage: H = FusionRing(G)
                sage: [x.q_dimension() for x in H.basis()]
                [1, 1, 1, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4]
                sage: sum(x.q_dimension()^2 for x in H.basis()) == G.order()^2
                True
            """
            if not self.is_simple_object():
                raise ValueError("Quantum dimension is only available for simple objects of a FusionRing")
            return self.parent().s_ij(self, self.parent().one(), base_coercion=base_coercion)

        def twist(self, reduced=True, base_coercion=True):
            r"""
            Return a rational number `h` such that `\theta = e^{i \pi h}`
            is the twist of ``self``. The quantity `e^{i \pi h}` is
            also available using :meth:`ribbon`.
            This method is only available for simple objects.

            EXAMPLES::

                sage: G21 = FusionRing("G2", 1)
                sage: [x.twist() for x in G21.basis()]
                [0, 4/5]
                sage: [G21.root_of_unity(x.twist()) for x in G21.basis()]
                [1, zeta60^14 - zeta60^4]
                sage: zeta60 = G21.field().gen()
                sage: zeta60^((4/5)*(60/2))
                zeta60^14 - zeta60^4

                sage: F42 = FusionRing("F4", 2)
                sage: [x.twist() for x in F42.basis()]
                [0, 18/11, 2/11, 12/11, 4/11]

                sage: E62 = FusionRing("E6", 2)
                sage: [x.twist() for x in E62.basis()]
                [0, 26/21, 12/7, 8/21, 8/21, 26/21, 2/3, 4/7, 2/3]

            EXAMPLES::

                sage: Q = FusionRing(CyclicPermutationGroup(3))
                sage: [x.twist() for x in Q.basis()]
                [0, 0, 0, 0, 2/3, 4/3, 0, 4/3, 2/3]
                sage: [x.ribbon() for x in Q.basis()]
                [1, 1, 1, 1, zeta3, -zeta3 - 1, 1, -zeta3 - 1, zeta3]
            """
            if not self.is_simple_object():
                raise ValueError("Quantum twist is only available for simple objects of a FusionRing")
            zeta = self.parent().field().gen()
            rib = self.ribbon()
            for k in range(4*self.parent()._cyclotomic_order):
                if zeta**k == rib:
                    return k/(2*self.parent()._cyclotomic_order)
