# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.graphs sage.modules
"""
Integrable Representations of Affine Lie Algebras
"""
# ***************************************************************************
#  Copyright (C) 2014, 2105 Daniel Bump <bump at match.stanford.edu>
#                           Travis Scrimshaw <tscrim at ucdavis.edu>
#                           Valentin Buciumas <buciumas at stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.category_object import CategoryObject
from sage.categories.modules import Modules
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import Matrix
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from sage.combinat.root_system.weyl_characters import WeylCharacterRing


# TODO: Make this a proper parent and implement actions
class IntegrableRepresentation(UniqueRepresentation, CategoryObject):
    r"""
    An irreducible integrable highest weight representation of
    an affine Lie algebra.

    INPUT:

    - ``Lam`` -- a dominant weight in an extended weight lattice
      of affine type

    REFERENCES:

    - [Ka1990]_

    .. [KMPS] Kass, Moody, Patera and Slansky, *Affine Lie algebras,
       weight multiplicities, and branching rules*. Vols. 1, 2. University
       of California Press, Berkeley, CA, 1990.

    .. [KacPeterson] Kac and Peterson. *Infinite-dimensional Lie algebras,
       theta functions and modular forms*. Adv. in Math. 53 (1984),
       no. 2, 125-264.

    .. [Carter] Carter, *Lie algebras of finite and affine type*. Cambridge
       University Press, 2005

    If `\Lambda` is a dominant integral weight for an affine root system,
    there exists a unique integrable representation `V=V_\Lambda` of highest
    weight `\Lambda`. If `\mu` is another weight, let `m(\mu)` denote the
    multiplicity of the weight `\mu` in this representation. The set
    `\operatorname{supp}(V)` of `\mu` such that `m(\mu) > 0` is contained in the
    paraboloid

    .. MATH::

         (\Lambda+\rho | \Lambda+\rho) - (\mu+\rho | \mu+\rho) \geq 0

    where `(\, | \,)` is the invariant inner product on the weight
    lattice and `\rho` is the Weyl vector. Moreover if `m(\mu)>0`
    then `\mu\in\operatorname{supp}(V)` differs from `\Lambda` by an element
    of the root lattice ([Ka1990]_, Propositions 11.3 and 11.4).

    Let `\delta` be the nullroot, which is the lowest positive imaginary
    root. Then by [Ka1990]_, Proposition 11.3 or Corollary 11.9, for fixed `\mu`
    the function `m(\mu - k\delta)` is a monotone increasing function of
    `k`. It is useful to take `\mu` to be such that this function is nonzero
    if and only if `k \geq 0`. Therefore we make the following definition.  If
    `\mu` is such that `m(\mu) \neq 0` but `m(\mu + \delta) = 0` then `\mu` is
    called *maximal*.

    Since `\delta` is fixed under the action of the affine Weyl group,
    and since the weight multiplicities are Weyl group invariant, the
    function `k \mapsto m(\mu - k \delta)` is unchanged if `\mu` is replaced
    by an equivalent weight. Therefore in tabulating these functions, we may
    assume that `\mu` is dominant. There are only a finite number of dominant
    maximal weights.

    Since every nonzero weight multiplicity appears in the string
    `\mu - k\delta` for one of the finite number of dominant maximal
    weights `\mu`, it is important to be able to compute these. We may
    do this as follows.

    EXAMPLES::

         sage: Lambda = RootSystem(['A',3,1]).weight_lattice(extended=true).fundamental_weights()
         sage: IntegrableRepresentation(Lambda[1]+Lambda[2]+Lambda[3]).print_strings()
         2*Lambda[0] + Lambda[2]: 4 31 161 665 2380 7658 22721 63120 166085 417295 1007601 2349655
         Lambda[0] + 2*Lambda[1]: 2 18 99 430 1593 5274 16005 45324 121200 308829 754884 1779570
         Lambda[0] + 2*Lambda[3]: 2 18 99 430 1593 5274 16005 45324 121200 308829 754884 1779570
         Lambda[1] + Lambda[2] + Lambda[3]: 1 10 60 274 1056 3601 11199 32354 88009 227555 563390 1343178
         3*Lambda[2] - delta: 3 21 107 450 1638 5367 16194 45687 121876 310056 757056 1783324
         sage: Lambda = RootSystem(['D',4,1]).weight_lattice(extended=true).fundamental_weights()
         sage: IntegrableRepresentation(Lambda[0]+Lambda[1]).print_strings()                        # long time
         Lambda[0] + Lambda[1]: 1 10 62 293 1165 4097 13120 38997 109036 289575 735870 1799620
         Lambda[3] + Lambda[4] - delta: 3 25 136 590 2205 7391 22780 65613 178660 463842 1155717 2777795

    In this example, we construct the extended weight lattice of Cartan
    type `A_3^{(1)}`, then define ``Lambda`` to be the fundamental
    weights `(\Lambda_i)_{i \in I}`. We find there are 5 maximal
    dominant weights in irreducible representation of highest weight
    `\Lambda_1 + \Lambda_2 + \Lambda_3`, and we determine their strings.

    It was shown in [KacPeterson]_ that each string is the set of Fourier
    coefficients of a modular form.

    Every weight `\mu` such that the weight multiplicity `m(\mu)` is
    nonzero has the form

      .. MATH::

          \Lambda - n_0 \alpha_0 - n_1 \alpha_1 - \cdots,

    where the `n_i` are nonnegative integers. This is represented internally
    as a tuple `(n_0, n_1, n_2, \ldots)`. If you want an individual
    multiplicity you use the method :meth:`m` and supply it with this tuple::

        sage: Lambda = RootSystem(['C',2,1]).weight_lattice(extended=true).fundamental_weights()
        sage: V = IntegrableRepresentation(2*Lambda[0]); V
        Integrable representation of ['C', 2, 1] with highest weight 2*Lambda[0]
        sage: V.m((3,5,3))
        18

    The :class:`IntegrableRepresentation` class has methods :meth:`to_weight`
    and :meth:`from_weight` to convert between this internal representation
    and the weight lattice::

        sage: delta = V.weight_lattice().null_root()
        sage: V.to_weight((4,3,2))
        -3*Lambda[0] + 6*Lambda[1] - Lambda[2] - 4*delta
        sage: V.from_weight(-3*Lambda[0] + 6*Lambda[1] - Lambda[2] - 4*delta)
        (4, 3, 2)

    To get more values, use the depth parameter::

        sage: L0 = RootSystem(["A",1,1]).weight_lattice(extended=true).fundamental_weight(0); L0
        Lambda[0]
        sage: IntegrableRepresentation(4*L0).print_strings(depth=20)
        4*Lambda[0]: 1 1 3 6 13 23 44 75 131 215 354 561 889 1368 2097 3153 4712 6936 10151 14677
        2*Lambda[0] + 2*Lambda[1] - delta: 1 2 5 10 20 36 66 112 190 310 501 788 1230 1880 2850 4256 6303 9222 13396 19262
        4*Lambda[1] - 2*delta: 1 2 6 11 23 41 75 126 215 347 561 878 1368 2082 3153 4690 6936 10121 14677 21055

    An example in type `C_2^{(1)}`::

        sage: Lambda = RootSystem(['C',2,1]).weight_lattice(extended=true).fundamental_weights()
        sage: V = IntegrableRepresentation(2*Lambda[0])
        sage: V.print_strings()    # long time
        2*Lambda[0]: 1 2 9 26 77 194 477 1084 2387 5010 10227 20198
        Lambda[0] + Lambda[2] - delta: 1 5 18 55 149 372 872 1941 4141 8523 17005 33019
        2*Lambda[1] - delta: 1 4 15 44 122 304 721 1612 3469 7176 14414 28124
        2*Lambda[2] - 2*delta: 2 7 26 72 194 467 1084 2367 5010 10191 20198 38907

    Examples for twisted affine types::

        sage: Lambda = RootSystem(["A",2,2]).weight_lattice(extended=True).fundamental_weights()
        sage: IntegrableRepresentation(Lambda[0]).strings()
        {Lambda[0]: [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56]}
        sage: Lambda = RootSystem(['G',2,1]).dual.weight_lattice(extended=true).fundamental_weights()
        sage: V = IntegrableRepresentation(Lambda[0]+Lambda[1]+Lambda[2])
        sage: V.print_strings() # long time
        6*Lambdacheck[0]: 4 28 100 320 944 2460 6064 14300 31968 69020 144676 293916
        3*Lambdacheck[0] + Lambdacheck[1]: 2 16 58 192 588 1568 3952 9520 21644 47456 100906 207536
        4*Lambdacheck[0] + Lambdacheck[2]: 4 22 84 276 800 2124 5288 12470 28116 61056 128304 261972
        2*Lambdacheck[1] - deltacheck: 2 8 32 120 354 980 2576 6244 14498 32480 69776 145528
        Lambdacheck[0] + Lambdacheck[1] + Lambdacheck[2]: 1 6 26 94 294 832 2184 5388 12634 28390 61488 128976
        2*Lambdacheck[0] + 2*Lambdacheck[2]: 2 12 48 164 492 1344 3428 8256 18960 41844 89208 184512
        3*Lambdacheck[2] - deltacheck: 4 16 60 208 592 1584 4032 9552 21728 47776 101068 207888
        sage: Lambda = RootSystem(['A',6,2]).weight_lattice(extended=true).fundamental_weights()
        sage: V = IntegrableRepresentation(Lambda[0]+2*Lambda[1])
        sage: V.print_strings() # long time
        5*Lambda[0]: 3 42 378 2508 13707 64650 272211 1045470 3721815 12425064 39254163 118191378
        3*Lambda[0] + Lambda[2]: 1 23 234 1690 9689 47313 204247 800029 2893198 9786257 31262198 95035357
        Lambda[0] + 2*Lambda[1]: 1 14 154 1160 6920 34756 153523 612354 2248318 7702198 24875351 76341630
        Lambda[0] + Lambda[1] + Lambda[3] - 2*delta: 6 87 751 4779 25060 113971 464842 1736620 6034717 19723537 61152367 181068152
        Lambda[0] + 2*Lambda[2] - 2*delta: 3 54 499 3349 18166 84836 353092 1341250 4725259 15625727 48938396 146190544
        Lambda[0] + 2*Lambda[3] - 4*delta: 15 195 1539 9186 45804 200073 789201 2866560 9723582 31120281 94724550 275919741
    """

    def __init__(self, Lam):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[1]+Lambda[2]+Lambda[3])

        Some methods required by the category are not implemented::

            sage: TestSuite(V).run()  # known bug (#21387)
        """
        CategoryObject.__init__(self, base=ZZ, category=Modules(ZZ))

        self._Lam = Lam
        self._P = Lam.parent()
        self._Q = self._P.root_system.root_lattice()

        # Store some extra simple computations that appear in tight loops
        self._Lam_rho = self._Lam + self._P.rho()

        self._cartan_matrix = self._P.root_system.cartan_matrix()
        self._cartan_type = self._P.root_system.cartan_type()

        self._classical_rank = self._cartan_type.classical().rank()
        self._index_set = self._P.index_set()
        self._index_set_classical = self._cartan_type.classical().index_set()
        self._cminv = self._cartan_type.classical().cartan_matrix().inverse()

        self._ddict = {}
        self._mdict = {tuple(0 for i in self._index_set): 1}
        # Coerce a classical root into the root lattice Q
        from_cl_root = lambda h: self._Q._from_dict(h._monomial_coefficients)
        self._classical_roots = [from_cl_root(al)
                                 for al in self._Q.classical().roots()]
        self._classical_positive_roots = [from_cl_root(al)
                                          for al in self._Q.classical().positive_roots()]
        self._a = self._cartan_type.a() # This is not cached
        self._ac = self._cartan_type.dual().a() # This is not cached
        self._eps = {i: self._a[i] / self._ac[i] for i in self._index_set}
        E = Matrix.diagonal([self._eps[i] for i in self._index_set_classical])
        self._ip = (self._cartan_type.classical().cartan_matrix()*E).inverse()

        # Extra data for the twisted cases
        if not self._cartan_type.is_untwisted_affine():
            self._classical_short_roots = frozenset(al for al in self._classical_roots
                                                    if self._inner_qq(al,al) == 2)

    def highest_weight(self):
        """
        Return the highest weight of ``self``.

        EXAMPLES::

            sage: Lambda = RootSystem(['D',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: IntegrableRepresentation(Lambda[0]+2*Lambda[2]).highest_weight()
            Lambda[0] + 2*Lambda[2]
        """
        return self._Lam

    def weight_lattice(self):
        """
        Return the weight lattice associated to ``self``.

        EXAMPLES::

            sage: V=IntegrableRepresentation(RootSystem(['E',6,1]).weight_lattice(extended=true).fundamental_weight(0))
            sage: V.weight_lattice()
            Extended weight lattice of the Root system of type ['E', 6, 1]
        """
        return self._P

    def root_lattice(self):
        """
        Return the root lattice associated to ``self``.

        EXAMPLES::

            sage: V=IntegrableRepresentation(RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weight(0))
            sage: V.root_lattice()
            Root lattice of the Root system of type ['F', 4, 1]
        """
        return self._Q

    @cached_method
    def level(self):
        r"""
        Return the level of ``self``.

        The level of a highest weight representation `V_{\Lambda}` is
        defined as `(\Lambda | \delta)` See [Ka1990]_ section 12.4.

        EXAMPLES::

            sage: Lambda = RootSystem(['G',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: [IntegrableRepresentation(Lambda[i]).level() for i in [0,1,2]]
            [1, 1, 2]
        """
        return ZZ(self._inner_pq(self._Lam, self._Q.null_root()))

    @cached_method
    def coxeter_number(self):
        """
        Return the Coxeter number of the Cartan type of ``self``.

        The Coxeter number is defined in [Ka1990]_ Chapter 6, and commonly
        denoted `h`.

        EXAMPLES::

            sage: Lambda = RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: V.coxeter_number()
            12
        """
        return sum(self._a)

    @cached_method
    def dual_coxeter_number(self):
        r"""
        Return the dual Coxeter number of the Cartan type of ``self``.

        The dual Coxeter number is defined in [Ka1990]_ Chapter 6, and commonly
        denoted `h^{\vee}`.

        EXAMPLES::

            sage: Lambda = RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: V.dual_coxeter_number()
            9
        """
        return sum(self._ac)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Lambda = RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: IntegrableRepresentation(Lambda[0])
            Integrable representation of ['F', 4, 1] with highest weight Lambda[0]
        """
        return "Integrable representation of %s with highest weight %s" % (self._cartan_type, self._Lam)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: Lambda = RootSystem(['C',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0]+2*Lambda[3])
            sage: latex(V)
            V_{\Lambda_{0} + 2 \Lambda_{3}}
        """
        return "V_{{{}}}".format(self._Lam._latex_())

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: Lambda = RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: V.cartan_type()
            ['F', 4, 1]
        """
        return self._cartan_type

    def _inner_qq(self, qelt1, qelt2):
        """
        Symmetric form between two elements of the root lattice
        associated to ``self``.

        EXAMPLES::

            sage: P = RootSystem(['F',4,1]).weight_lattice(extended=true)
            sage: Lambda = P.fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: alpha = V.root_lattice().simple_roots()
            sage: Matrix([[V._inner_qq(alpha[i], alpha[j]) for j in V._index_set] for i in V._index_set])
            [   2   -1    0    0    0]
            [  -1    2   -1    0    0]
            [   0   -1    2   -1    0]
            [   0    0   -1    1 -1/2]
            [   0    0    0 -1/2    1]

        .. WARNING::

            If ``qelt1`` or ``qelt1`` accidentally gets coerced into
            the extended weight lattice, this will return an answer,
            and it will be wrong. To make this code robust, parents
            should be checked. This is not done since in the application
            the parents are known, so checking would unnecessarily slow
            us down.
        """
        mc1 = qelt1.monomial_coefficients()
        mc2 = qelt2.monomial_coefficients()
        zero = ZZ.zero()
        return sum(mc1.get(i, zero) * mc2.get(j, zero)
                   * self._cartan_matrix[i,j] / self._eps[i]
                   for i in self._index_set for j in self._index_set)

    def _inner_pq(self, pelt, qelt):
        """
        Symmetric form between an element of the weight and root lattices
        associated to ``self``.

        .. WARNING::

            If ``qelt`` accidentally gets coerced into the extended weight
            lattice, this will return an answer, and it will be wrong. To make
            this code robust, parents should be checked. This is not done
            since in the application the parents are known, so checking would
            unnecessarily slow us down.

        EXAMPLES::

            sage: P = RootSystem(['F',4,1]).weight_lattice(extended=true)
            sage: Lambda = P.fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: alpha = V.root_lattice().simple_roots()
            sage: Matrix([[V._inner_pq(P(alpha[i]), alpha[j]) for j in V._index_set] for i in V._index_set])
            [   2   -1    0    0    0]
            [  -1    2   -1    0    0]
            [   0   -1    2   -1    0]
            [   0    0   -1    1 -1/2]
            [   0    0    0 -1/2    1]
            sage: P = RootSystem(['G',2,1]).weight_lattice(extended=true)
            sage: P = RootSystem(['G',2,1]).weight_lattice(extended=true)
            sage: Lambda = P.fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: alpha = V.root_lattice().simple_roots()
            sage: Matrix([[V._inner_pq(Lambda[i],alpha[j]) for j in V._index_set] for i in V._index_set])
            [  1   0   0]
            [  0 1/3   0]
            [  0   0   1]
        """
        mcp = pelt.monomial_coefficients()
        mcq = qelt.monomial_coefficients()
        zero = ZZ.zero()
        return sum(mcp.get(i, zero) * mcq[i] / self._eps[i] for i in mcq)

    def _inner_pp(self, pelt1, pelt2):
        """
        Symmetric form between two elements of the weight lattice
        associated to ``self``.

        EXAMPLES::

            sage: P = RootSystem(['G',2,1]).weight_lattice(extended=true)
            sage: Lambda = P.fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: alpha = V.root_lattice().simple_roots()
            sage: Matrix([[V._inner_pp(Lambda[i],P(alpha[j])) for j in V._index_set] for i in V._index_set])
            [  1   0   0]
            [  0 1/3   0]
            [  0   0   1]
            sage: Matrix([[V._inner_pp(Lambda[i],Lambda[j]) for j in V._index_set] for i in V._index_set])
            [  0   0   0]
            [  0 2/3   1]
            [  0   1   2]
        """
        mc1 = pelt1.monomial_coefficients()
        mc2 = pelt2.monomial_coefficients()
        zero = ZZ.zero()
        mc1d = mc1.get('delta', zero)
        mc2d = mc2.get('delta', zero)
        return sum(mc1.get(i,zero) * self._ac[i] * mc2d
                   + mc2.get(i,zero) * self._ac[i] * mc1d
                   for i in self._index_set) \
               + sum(mc1.get(i,zero) * mc2.get(j,zero) * self._ip[ii,ij]
                     for ii, i in enumerate(self._index_set_classical)
                     for ij, j in enumerate(self._index_set_classical))

    def to_weight(self, n):
        r"""
        Return the weight associated to the tuple ``n`` in ``self``.

        If ``n`` is the tuple `(n_1, n_2, \ldots)`, then the associated
        weight is `\Lambda - \sum_i n_i \alpha_i`, where `\Lambda`
        is the weight of the representation.

        INPUT:

        - ``n`` -- tuple representing a weight

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(2*Lambda[2])
            sage: V.to_weight((1,0,0))
            -2*Lambda[0] + Lambda[1] + 3*Lambda[2] - delta
        """
        alpha = self._P.simple_roots()
        I = self._index_set
        return self._Lam - self._P.sum(val * alpha[I[i]]
                                       for i,val in enumerate(n))

    def _from_weight_helper(self, mu, check=False):
        r"""
        Return the coefficients of a tuple of the weight ``mu`` expressed
        in terms of the simple roots in ``self``.

        The tuple ``n`` is defined as the tuple `(n_0, n_1, \ldots)`
        such that `\mu = \sum_{i \in I} n_i \alpha_i`.

        INPUT:

        - ``mu`` -- an element in the root lattice

        .. TODO::

            Implement this as a section map of the inverse of the
            coercion from `Q \to P`.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(2*Lambda[2])
            sage: V.to_weight((1,0,0))
            -2*Lambda[0] + Lambda[1] + 3*Lambda[2] - delta
            sage: delta = V.weight_lattice().null_root()
            sage: V._from_weight_helper(2*Lambda[0] - Lambda[1] - 1*Lambda[2] + delta)
            (1, 0, 0)
        """
        mu = self._P(mu)
        zero = ZZ.zero()
        n0 = mu.monomial_coefficients().get('delta', zero)
        mu0 = mu - n0 * self._P.simple_root(self._cartan_type.special_node())
        ret = [n0] # This should be in ZZ because it is in the weight lattice
        mc_mu0 = mu0.monomial_coefficients()
        for ii, i in enumerate(self._index_set_classical):
            # -1 for indexing
            ret.append( sum(self._cminv[ii,ij] * mc_mu0.get(j, zero)
                               for ij, j in enumerate(self._index_set_classical)) )
        if check:
            return all(x in ZZ for x in ret)
        else:
            return tuple(ZZ(x) for x in ret)

    def from_weight(self, mu):
        r"""
        Return the tuple `(n_0, n_1, ...)`` such that ``mu`` equals
        `\Lambda - \sum_{i \in I} n_i \alpha_i` in ``self``, where `\Lambda`
        is the highest weight of ``self``.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(2*Lambda[2])
            sage: V.to_weight((1,0,0))
            -2*Lambda[0] + Lambda[1] + 3*Lambda[2] - delta
            sage: delta = V.weight_lattice().null_root()
            sage: V.from_weight(-2*Lambda[0] + Lambda[1] + 3*Lambda[2] - delta)
            (1, 0, 0)
        """
        return self._from_weight_helper(self._Lam - mu)

    def s(self, n, i):
        """
        Return the action of the ``i``-th simple reflection on the
        internal representation of weights by tuples ``n`` in ``self``.

        EXAMPLES::

            sage: V = IntegrableRepresentation(RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weight(0))
            sage: [V.s((0,0,0),i) for i in V._index_set]
            [(1, 0, 0), (0, 0, 0), (0, 0, 0)]
        """
        ret = list(n) # This makes a copy
        ret[i] += self._Lam._monomial_coefficients.get(i, ZZ.zero())
        ret[i] -= sum(val * self._cartan_matrix[i,j] for j,val in enumerate(n))
        return tuple(ret)

    def to_dominant(self, n):
        """
        Return the dominant weight in ``self`` equivalent to ``n``
        under the affine Weyl group.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(3*Lambda[0])
            sage: n = V.to_dominant((13,11,7)); n
            (4, 3, 3)
            sage: V.to_weight(n)
            Lambda[0] + Lambda[1] + Lambda[2] - 4*delta
        """
        if n in self._ddict:
            return self._ddict[n]

        path = [n]
        alpha = self._P.simple_roots()
        next = True
        cur_wt = self.to_weight(n)

        while next:
            if path[-1] in self._ddict:
                path.append( self._ddict[path[-1]] )
                break

            next = False
            mc = cur_wt.monomial_coefficients()
            # Most weights are dense over the index set
            for i in self._index_set:
                if mc.get(i, 0) < 0:
                    m = self.s(path[-1], i)
                    if m in self._ddict:
                        path.append(self._ddict[m])
                    else:
                        cur_wt -= (m[i] - path[-1][i]) * alpha[i]
                        path.append(m)
                        next = True
                    break

        # We don't want any dominant weight to refer to itself in self._ddict
        #   as this leads to an infinite loop with self.m() when the dominant
        #   weight does not have a known multiplicity.
        v = path.pop()
        for m in path:
            self._ddict[m] = v
        return v

    def _freudenthal_roots_imaginary(self, nu):
        r"""
        Iterate over the set of imaginary roots `\alpha \in \Delta^+`
        in ``self`` such that `\nu - \alpha \in Q^+`.

        INPUT:

        - ``nu`` -- an element in `Q`

        EXAMPLES::

            sage: Lambda = RootSystem(['B',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0]+Lambda[1]+Lambda[3])
            sage: [list(V._freudenthal_roots_imaginary(V.highest_weight() - mw))
            ....:  for mw in V.dominant_maximal_weights()]
            [[], [], [], [], []]
        """
        l = self._from_weight_helper(nu)
        kp = min(l[i] // self._a[i] for i in self._index_set)
        delta = self._Q.null_root()
        for u in range(1, kp+1):
            yield u * delta

    def _freudenthal_roots_real(self, nu):
        r"""
        Iterate over the set of real positive roots `\alpha \in \Delta^+`
        in ``self`` such that `\nu - \alpha \in Q^+`.

        See [Ka1990]_ Proposition 6.3 for the way to compute the set of real
        roots for twisted affine case.

        INPUT:

        - ``nu`` -- an element in `Q`

        EXAMPLES::

            sage: Lambda = RootSystem(['B',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0]+Lambda[1]+Lambda[3])
            sage: mw = V.dominant_maximal_weights()[0]
            sage: sorted(V._freudenthal_roots_real(V.highest_weight() - mw), key=str)
            [alpha[1],
             alpha[1] + alpha[2],
             alpha[1] + alpha[2] + alpha[3],
             alpha[2],
             alpha[2] + alpha[3],
             alpha[3]]
        """
        for al in self._classical_positive_roots:
            if min(self._from_weight_helper(nu-al)) >= 0:
                yield al

        if self._cartan_type.is_untwisted_affine():
            # untwisted case
            for al in self._classical_roots:
                for ir in self._freudenthal_roots_imaginary(nu-al):
                    yield al + ir

        elif self._cartan_type.type() == 'BC':
            #case A^2_{2l}
            # We have to keep track of the roots we have visited for this case
            ret = set(self._classical_positive_roots)
            for al in self._classical_roots:
                if al in self._classical_short_roots:
                    for ir in self._freudenthal_roots_imaginary(nu-al):
                        ret.add(al + ir)
                        yield al + ir
                else:
                    fri = list(self._freudenthal_roots_imaginary(nu-al))
                    friset = set(fri)
                    for ir in fri:
                        if 2*ir in friset:
                            ret.add(al + 2*ir)
                            yield al + 2*ir
                    alpha = self._Q.simple_roots()
                    fri = list(self._freudenthal_roots_imaginary(2*nu-al))
                    for ir in fri[::2]:
                        rt = sum( val // 2 * alpha[i] for i,val in
                                  enumerate(self._from_weight_helper(al+ir)) )
                        if rt not in ret:
                            ret.add(rt)
                            yield rt

        elif self._cartan_type.dual().type() == 'G':
            # case D^3_4 in the Kac notation
            for al in self._classical_roots:
                if al in self._classical_short_roots:
                    for ir in self._freudenthal_roots_imaginary(nu-al):
                        yield al + ir
                else:
                    fri = list(self._freudenthal_roots_imaginary(nu-al))
                    friset = set(fri)
                    for ir in fri:
                        if 3*ir in friset:
                            yield al + 3*ir

        elif self._cartan_type.dual().type() in ['B','C','F']:
            #case A^2_{2l-1} or case D^2_{l+1} or case E^2_6:
            for al in self._classical_roots:
                if al in self._classical_short_roots:
                    for ir in self._freudenthal_roots_imaginary(nu-al):
                        yield al + ir
                else:
                    fri = list(self._freudenthal_roots_imaginary(nu-al))
                    friset = set(fri)
                    for ir in fri:
                        if 2*ir in friset:
                            yield al + 2*ir

    def _freudenthal_accum(self, nu, al):
        """
        Helper method for computing the Freudenthal formula in ``self``.

        EXAMPLES::

            sage: Lambda = RootSystem(['B',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0]+Lambda[1]+Lambda[3])
            sage: mw = V.dominant_maximal_weights()[0]
            sage: F = V._freudenthal_roots_real(V.highest_weight() - mw)
            sage: sorted([V._freudenthal_accum(mw, al) for al in F])
            [3, 3, 3, 4, 4, 4]
        """
        ret = 0
        n = list(self._from_weight_helper(self._Lam - nu))
        ip = self._inner_pq(nu, al)
        n_shift = self._from_weight_helper(al)
        ip_shift = self._inner_qq(al, al)

        while min(n) >= 0:
            # Change in data by adding ``al`` to our current weight
            ip += ip_shift
            for i,val in enumerate(n_shift):
                n[i] -= val
            # Compute the multiplicity
            ret += 2 * self.m(tuple(n)) * ip
        return ret

    def _m_freudenthal(self, n):
        r"""
        Compute the weight multiplicity using the Freudenthal
        multiplicity formula in ``self``.

        The multiplicities of the imaginary roots for the twisted
        affine case are different than those for the untwisted case.
        See [Carter]_ Corollary 18.10 for general type and Corollary
        18.15 for `A^2_{2l}`

        EXAMPLES::

            sage: Lambda = RootSystem(['B',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0]+Lambda[1]+Lambda[3])
            sage: D = list(V.dominant_maximal_weights())
            sage: D.remove(V.highest_weight())
            sage: [V._m_freudenthal(V.from_weight(mw)) for mw in D]
            [3, 7, 3, 3]
        """
        if min(n) < 0:
            return 0
        mu = self.to_weight(n)
        I = self._index_set
        al = self._Q._from_dict({I[i]: val for i,val in enumerate(n) if val},
                                remove_zeros=False)
        cr = self._classical_rank
        num = sum(self._freudenthal_accum(mu, alr)
                  for alr in self._freudenthal_roots_real(self._Lam - mu))

        if self._cartan_type.is_untwisted_affine():
            num += sum(cr * self._freudenthal_accum(mu, alr)
                       for alr in self._freudenthal_roots_imaginary(self._Lam - mu))

        elif self._cartan_type.dual().type() == 'B': # A_{2n-1}^{(2)}
            val = 1
            for rt in self._freudenthal_roots_imaginary(self._Lam - mu):
                # k-th element (starting from 1) is k*delta
                num += (cr - val) * self._freudenthal_accum(mu, rt)
                val = 1 - val

        elif self._cartan_type.type() == 'BC': # A_{2n}^{(2)}
            num += sum(cr * self._freudenthal_accum(mu, alr)
                       for alr in self._freudenthal_roots_imaginary(self._Lam - mu))

        elif self._cartan_type.dual() == 'C': # D_{n+1}^{(2)}
            val = 1
            for rt in self._freudenthal_roots_imaginary(self._Lam - mu):
                # k-th element (starting from 1) is k*delta
                num += (cr - (cr - 1)*val) * self._freudenthal_accum(mu, rt)
                val = 1 - val

        elif self._cartan_type.dual().type() == 'F': # E_6^{(2)}
            val = 1
            for rt in self._freudenthal_roots_imaginary(self._Lam - mu):
                # k-th element (starting from 1) is k*delta
                num += (4 - 2*val) * self._freudenthal_accum(mu, rt)
                val = 1 - val

        elif self._cartan_type.dual().type() == 'G': # D_4^{(3)} (or dual of G_2^{(1)})
            for k,rt in enumerate(self._freudenthal_roots_imaginary(self._Lam - mu)):
                # k-th element (starting from 1) is k*delta
                if (k+1) % 3 == 0:
                    num += 2 * self._freudenthal_accum(mu, rt)
                else:
                    num += self._freudenthal_accum(mu, rt)

        den = 2*self._inner_pq(self._Lam_rho, al) - self._inner_qq(al, al)
        try:
            return ZZ(num / den)
        except TypeError:
            return None

    def m(self, n):
        r"""
        Return the multiplicity of the weight `\mu` in ``self``, where
        `\mu = \Lambda - \sum_i n_i \alpha_i`.

        INPUT:

        - ``n`` -- tuple representing a weight `\mu`

        EXAMPLES::

            sage: Lambda = RootSystem(['E',6,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: u = V.highest_weight() - V.weight_lattice().null_root()
            sage: V.from_weight(u)
            (1, 1, 2, 2, 3, 2, 1)
            sage: V.m(V.from_weight(u))
            6
        """
        # TODO: Make this non-recursive by implementing our own stack
        # The recursion follows:
        #   - m
        #   - _m_freudenthal
        #   - _freudenthal_accum
        if n in self._mdict:
            return self._mdict[n]
        elif n in self._ddict:
            self._mdict[n] = self.m(self._ddict[n])
        m = self.to_dominant(n)
        if m in self._mdict:
            return self._mdict[m]
        ret = self._m_freudenthal(m)
        assert ret is not None, "m: error - failed to compute m{}".format(n)
        self._mdict[n] = ret
        return ret

    def mult(self, mu):
        """
        Return the weight multiplicity of ``mu``.

        INPUT:

        - ``mu`` -- an element of the weight lattice

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: L = RootSystem("B3~").weight_lattice(extended=True)
            sage: Lambda = L.fundamental_weights()
            sage: delta = L.null_root()
            sage: W = L.weyl_group(prefix='s')
            sage: [s0,s1,s2,s3] = W.simple_reflections()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: V.mult(Lambda[2] - 2*delta)
            3
            sage: V.mult(Lambda[2] - Lambda[1])
            0
            sage: weights = [w.action(Lambda[1] - 4*delta) for w in [s1,s2,s0*s1*s2*s3]]
            sage: weights
            [-Lambda[1] + Lambda[2] - 4*delta,
             Lambda[1] - 4*delta,
             -Lambda[1] + Lambda[2] - 4*delta]
            sage: [V.mult(mu) for mu in weights]
            [35, 35, 35]

        TESTS::

            sage: L = RootSystem("B3~").weight_lattice(extended=True)
            sage: La = L.fundamental_weights()
            sage: V = IntegrableRepresentation(La[0])
            sage: Q = RootSystem("B3~").root_space()
            sage: al = Q.simple_roots()
            sage: V.mult(1/2*al[1])
            0
        """
        try:
            n = self.from_weight(mu)
        except TypeError:
            return ZZ.zero()
        return self.m(n)

    @cached_method
    def dominant_maximal_weights(self):
        r"""
        Return the dominant maximal weights of ``self``.

        A weight `\mu` is *maximal* if it has nonzero multiplicity but
        `\mu + \delta`` has multiplicity zero. There are a finite number
        of dominant maximal weights. Indeed, [Ka1990]_ Proposition 12.6
        shows that the dominant maximal weights are in bijection with
        the classical weights in `k \cdot F` where `F` is the fundamental
        alcove and `k` is the level. The construction used in this
        method is based on that Proposition.

        EXAMPLES::

            sage: Lambda = RootSystem(['C',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: IntegrableRepresentation(2*Lambda[0]).dominant_maximal_weights()
            (2*Lambda[0],
             Lambda[0] + Lambda[2] - delta,
             2*Lambda[1] - delta,
             Lambda[1] + Lambda[3] - 2*delta,
             2*Lambda[2] - 2*delta,
             2*Lambda[3] - 3*delta)
        """
        k = self.level()
        Lambda = self._P.fundamental_weights()

        def next_level(wt):
            return [wt + Lambda[i] for i in self._index_set_classical
                    if (wt + Lambda[i]).level() <= k]
        R = RecursivelyEnumeratedSet([self._P.zero()], next_level)
        candidates = [x + (k - x.level())*Lambda[0] for x in list(R)]
        ret = []
        delta = self._Q.null_root()
        for x in candidates:
            if self._from_weight_helper(self._Lam-x, check=True):
                t = 0
                while self.m(self.from_weight(x - t*delta)) == 0:
                    t += 1
                ret.append(x - t*delta)
        return tuple(ret)

    def string(self, max_weight, depth=12):
        r"""
        Return the list of multiplicities `m(\Lambda - k \delta)` in
        ``self``, where `\Lambda` is ``max_weight`` and `k` runs from `0`
        to ``depth``.

        INPUT:

        - ``max_weight`` -- a dominant maximal weight
        - ``depth`` -- (default: 12) the maximum value of `k`

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(2*Lambda[0])
            sage: V.string(2*Lambda[0])
            [1, 2, 8, 20, 52, 116, 256, 522, 1045, 1996, 3736, 6780]
            sage: V.string(Lambda[1] + Lambda[2])
            [0, 1, 4, 12, 32, 77, 172, 365, 740, 1445, 2736, 5041]
        """
        ret = []
        delta = self._Q.null_root()
        cur_weight = max_weight
        for k in range(depth):
            ret.append(self.m( self.from_weight(cur_weight) ))
            cur_weight -= delta
        return ret

    def strings(self, depth=12):
        """
        Return the set of dominant maximal weights of ``self``, together
        with the string coefficients for each.

        OPTIONAL:

        - ``depth`` -- (default: 12) a parameter indicating how far
          to push computations

        EXAMPLES::

            sage: Lambda = RootSystem(['A',1,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(2*Lambda[0])
            sage: S = V.strings(depth=25)
            sage: for k in S:
            ....:     print("{}: {}".format(k, ' '.join(str(x) for x in S[k])))
            2*Lambda[0]: 1 1 3 5 10 16 28 43 70 105 161 236 350 501 722 1016 1431 1981 2741 3740 5096 6868 9233 12306 16357
            2*Lambda[1] - delta: 1 2 4 7 13 21 35 55 86 130 196 287 420 602 858 1206 1687 2331 3206 4368 5922 7967 10670 14193 18803
        """
        return {max_weight: self.string(max_weight, depth)
                for max_weight in self.dominant_maximal_weights()}

    def print_strings(self, depth=12):
        """
        Print the strings of ``self``.

        .. SEEALSO::

            :meth:`strings`

        EXAMPLES::

            sage: Lambda = RootSystem(['A',1,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(2*Lambda[0])
            sage: V.print_strings(depth=25)
            2*Lambda[0]: 1 1 3 5 10 16 28 43 70 105 161 236 350 501 722 1016 1431 1981 2741 3740 5096 6868 9233 12306 16357
            2*Lambda[1] - delta: 1 2 4 7 13 21 35 55 86 130 196 287 420 602 858 1206 1687 2331 3206 4368 5922 7967 10670 14193 18803
        """
        S = self.strings(depth=depth)
        for mw in self.dominant_maximal_weights():
            print("{}: {}".format(mw, ' '.join(str(x) for x in S[mw])) )

    def modular_characteristic(self, mu=None):
        r"""
        Return the modular characteristic of ``self``.

        The modular characteristic is a rational number introduced
        by Kac and Peterson [KacPeterson]_, required to interpret the
        string functions as Fourier coefficients of modular forms. See
        [Ka1990]_ Section 12.7. Let `k` be the level, and let `h^\vee`
        be the dual Coxeter number. Then

        .. MATH::

            m_\Lambda = \frac{|\Lambda+\rho|^2}{2(k+h^\vee)}
            - \frac{|\rho|^2}{2h^\vee}

        If `\mu` is a weight, then

        .. MATH::

            m_{\Lambda,\mu} = m_\Lambda - \frac{|\mu|^2}{2k}.

        OPTIONAL:

        - ``mu`` -- a weight, or alternatively,
        - ``n`` -- tuple representing a weight `\mu`

        If no optional parameter is specified, this returns `m_\Lambda`.
        If ``mu`` is specified, it returns `m_{\Lambda,\mu}`. You may
        use the tuple ``n`` to specify `\mu`. If you do this, `\mu` is
        `\Lambda - \sum_i n_i \alpha_i`.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',1,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(3*Lambda[0]+2*Lambda[1])
            sage: [V.modular_characteristic(x) for x in V.dominant_maximal_weights()]
            [11/56, -1/280, 111/280]
        """
        if isinstance(mu, tuple):
            n = mu
        else:
            n = self.from_weight(mu)
        k = self.level()
        hd = self.dual_coxeter_number()
        rho = self._P.rho()
        m_Lambda = self._inner_pp(self._Lam_rho, self._Lam_rho) / (2*(k+hd)) \
                   - self._inner_pp(rho, rho) / (2*hd)
        if n is None:
            return m_Lambda
        mu = self.to_weight(n)
        return m_Lambda - self._inner_pp(mu,mu) / (2*k)

    def branch(self, i=None, weyl_character_ring=None, sequence=None, depth=5):
        r"""
        Return the branching rule on ``self``.

        Removing any node from the extended Dynkin diagram of the affine
        Lie algebra results in the Dynkin diagram of a classical Lie
        algebra, which is therefore a Lie subalgebra. For example
        removing the `0` node from the Dynkin diagram of type ``[X, r, 1]``
        produces the classical Dynkin diagram of ``[X, r]``.

        Thus for each `i` in the index set, we may restrict ``self`` to
        the corresponding classical subalgebra. Of course ``self`` is
        an infinite dimensional representation, but each weight `\mu`
        is assigned a grading by the number of times the simple root
        `\alpha_i` appears in `\Lambda-\mu`. Thus the branched
        representation is graded and we get sequence of finite-dimensional
        representations which this method is able to compute.

        OPTIONAL:

        - ``i`` -- (default: 0) an element of the index set
        - ``weyl_character_ring`` -- a WeylCharacterRing
        - ``sequence`` -- dictionary
        - ``depth`` -- (default: 5) an upper bound for `k` determining
          how many terms to give

        In the default case where `i = 0`, you do not need to specify
        anything else, though you may want to increase the depth if
        you need more terms.

        EXAMPLES::

            sage: Lambda = RootSystem(['A',2,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(2*Lambda[0])
            sage: b = V.branch(); b                                                     # needs sage.libs.gap
            [A2(0,0),
             A2(1,1),
             A2(0,0) + 2*A2(1,1) + A2(2,2),
             2*A2(0,0) + 2*A2(0,3) + 4*A2(1,1) + 2*A2(3,0) + 2*A2(2,2),
             4*A2(0,0) + 3*A2(0,3) + 10*A2(1,1) + 3*A2(3,0) + A2(1,4) + 6*A2(2,2) + A2(4,1),
             6*A2(0,0) + 9*A2(0,3) + 20*A2(1,1) + 9*A2(3,0) + 3*A2(1,4) + 12*A2(2,2) + 3*A2(4,1) + A2(3,3)]

        If the parameter ``weyl_character_ring`` is omitted, the ring may be recovered
        as the parent of one of the branched coefficients::

            sage: A2 = b[0].parent(); A2                                                # needs sage.libs.gap
            The Weyl Character Ring of Type A2 with Integer Ring coefficients

        If `i` is not zero then you should specify the :class:`WeylCharacterRing` that you
        are branching to. This is determined by the Dynkin diagram::

            sage: Lambda = RootSystem(['B',3,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: V.cartan_type().dynkin_diagram()
                O 0
                |
                |
            O---O=>=O
            1   2   3
            B3~

        In this example, we observe that removing the `i=2` node from the
        Dynkin diagram produces a reducible diagram of type ``A1xA1xA1``.
        Thus we have a branching to
        `\mathfrak{sl}(2) \times \mathfrak{sl}(2) \times \mathfrak{sl}(2)`::

            sage: A1xA1xA1 = WeylCharacterRing("A1xA1xA1",style='coroots')              # needs sage.libs.gap
            sage: V.branch(i=2,weyl_character_ring=A1xA1xA1)                            # needs sage.libs.gap
            [A1xA1xA1(1,0,0),
             A1xA1xA1(0,1,2),
             A1xA1xA1(1,0,0) + A1xA1xA1(1,2,0) + A1xA1xA1(1,0,2),
             A1xA1xA1(2,1,2) + A1xA1xA1(0,1,0) + 2*A1xA1xA1(0,1,2),
             3*A1xA1xA1(1,0,0) + 2*A1xA1xA1(1,2,0) + A1xA1xA1(1,2,2) + 2*A1xA1xA1(1,0,2) + A1xA1xA1(1,0,4) + A1xA1xA1(3,0,0),
             A1xA1xA1(2,1,0) + 3*A1xA1xA1(2,1,2) + 2*A1xA1xA1(0,1,0) + 5*A1xA1xA1(0,1,2) + A1xA1xA1(0,1,4) + A1xA1xA1(0,3,2)]

        If the nodes of the two Dynkin diagrams are not in the same order, you
        must specify an additional parameter, ``sequence`` which gives a dictionary
        to the affine Dynkin diagram to the classical one.

        EXAMPLES::

            sage: Lambda = RootSystem(['F',4,1]).weight_lattice(extended=true).fundamental_weights()
            sage: V = IntegrableRepresentation(Lambda[0])
            sage: V.cartan_type().dynkin_diagram()
            O---O---O=>=O---O
            0   1   2   3   4
            F4~
            sage: A1xC3=WeylCharacterRing("A1xC3",style='coroots')
            sage: A1xC3.dynkin_diagram()
            O
            1
            O---O=<=O
            2   3   4
            A1xC3

        Observe that removing the `i=1` node from the ``F4~`` Dynkin diagram
        gives the ``A1xC3`` diagram, but the roots are in a different order.
        The nodes `0, 2, 3, 4` of ``F4~`` correspond to ``1, 4, 3, 2``
        of ``A1xC3`` and so we encode this in a dictionary::

            sage: V.branch(i=1, weyl_character_ring=A1xC3, sequence={0:1,2:4,3:3,4:2})  # long time
            [A1xC3(1,0,0,0),
             A1xC3(0,0,0,1),
             A1xC3(1,0,0,0) + A1xC3(1,2,0,0),
             A1xC3(2,0,0,1) + A1xC3(0,0,0,1) + A1xC3(0,1,1,0),
             2*A1xC3(1,0,0,0) + A1xC3(1,0,1,0) + 2*A1xC3(1,2,0,0) + A1xC3(1,0,2,0) + A1xC3(3,0,0,0),
             2*A1xC3(2,0,0,1) + A1xC3(2,1,1,0) + A1xC3(0,1,0,0) + 3*A1xC3(0,0,0,1) + 2*A1xC3(0,1,1,0) + A1xC3(0,2,0,1)]

        The branch method gives a way of computing the graded dimension of the integrable representation::

            sage: Lambda = RootSystem("A1~").weight_lattice(extended=true).fundamental_weights()
            sage: V=IntegrableRepresentation(Lambda[0])
            sage: r = [x.degree() for x in V.branch(depth=15)]; r
            [1, 3, 4, 7, 13, 19, 29, 43, 62, 90, 126, 174, 239, 325, 435, 580]
            sage: oeis(r)                                                        # optional -- internet
            0: A029552: Expansion of phi(x) / f(-x) in powers of x where phi(), f() are Ramanujan theta functions.
        """
        if i is None:
            i = self._cartan_type.special_node()
        if i == self._cartan_type.special_node() or self._cartan_type.type() == 'A':
            if weyl_character_ring is None:
                weyl_character_ring = WeylCharacterRing(self._cartan_type.classical(), style='coroots')
            if weyl_character_ring.cartan_type() != self._cartan_type.classical():
                raise ValueError("Cartan type of WeylCharacterRing must be %s" % self.cartan_type().classical())
        elif weyl_character_ring is None:
            raise ValueError("the argument weyl_character_ring cannot be omitted if i != 0")
        if sequence is None:
            sequence = {}
            for j in self._index_set:
                if j < i:
                    sequence[j] = j+1
                elif j > i:
                    sequence[j] = j

        def next_level(x):
            ret = []
            for j in self._index_set:
                t = list(x[0])
                t[j] += 1
                t = tuple(t)
                m = self.m(t)
                if m > 0 and t[i] <= depth:
                    ret.append((t,m))
            return ret
        hwv = (tuple([0 for j in self._index_set]), 1)
        terms = RecursivelyEnumeratedSet([hwv], next_level)
        fw = weyl_character_ring.fundamental_weights()
        P = self.weight_lattice()
        ret = []
        for l in range(depth+1):
            lterms = [x for x in terms if x[0][i] == l]
            ldict = {}
            for x in lterms:
                mc = P(self.to_weight(x[0])).monomial_coefficients()
                contr = sum(fw[sequence[j]]*mc.get(j,0)
                            for j in self._index_set if j != i).coerce_to_sl()
                if contr in ldict:
                    ldict[contr] += x[1]
                else:
                    ldict[contr] = x[1]
            ret.append(weyl_character_ring.char_from_weights(ldict))
        return ret
