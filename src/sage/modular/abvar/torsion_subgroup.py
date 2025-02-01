# sage.doctest: needs sage.libs.flint sage.libs.pari
"""
Torsion subgroups of modular abelian varieties

Sage can compute information about the structure of the torsion
subgroup of a modular abelian variety. Sage computes a multiple of
the order by computing the greatest common divisor of the orders of
the torsion subgroup of the reduction of the abelian variety modulo
p for various primes p. Sage computes a divisor of the order by
computing the rational cuspidal subgroup. When these two bounds
agree (which is often the case), we determine the exact structure
of the torsion subgroup.

AUTHORS:

- William Stein (2007-03)

EXAMPLES: First we consider `J_0(50)` where everything
works out nicely::

    sage: J = J0(50)
    sage: T = J.rational_torsion_subgroup(); T
    Torsion subgroup of Abelian variety J0(50) of dimension 2
    sage: T.multiple_of_order()
    15
    sage: T.divisor_of_order()
    15
    sage: T.gens()
    [[(1/15, 3/5, 2/5, 14/15)]]
    sage: T.invariants()
    [15]
    sage: d = J.decomposition(); d
    [Simple abelian subvariety 50a(1,50) of dimension 1 of J0(50),
     Simple abelian subvariety 50b(1,50) of dimension 1 of J0(50)]
    sage: d[0].rational_torsion_subgroup().order()
    3
    sage: d[1].rational_torsion_subgroup().order()
    5

Next we make a table of the upper and lower bounds for each new
factor.

::

    sage: for N in range(1,38):
    ....:    for A in J0(N).new_subvariety().decomposition():
    ....:        T = A.rational_torsion_subgroup()
    ....:        print('%-5s%-5s%-5s%-5s'%(N, A.dimension(), T.divisor_of_order(), T.multiple_of_order()))
    11   1    5    5
    14   1    6    6
    15   1    8    8
    17   1    4    4
    19   1    3    3
    20   1    6    6
    21   1    8    8
    23   2    11   11
    24   1    8    8
    26   1    3    3
    26   1    7    7
    27   1    3    3
    29   2    7    7
    30   1    6    6
    31   2    5    5
    32   1    4    4
    33   1    4    4
    34   1    6    6
    35   1    3    3
    35   2    16   16
    36   1    6    6
    37   1    1    1
    37   1    3    3

TESTS::

    sage: T = J0(54).rational_torsion_subgroup()
    sage: loads(dumps(T)) == T
    True
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import divisors, gcd
from sage.misc.misc_c import prod
from sage.modular.abvar.torsion_point import TorsionPoint
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_class
from sage.modular.arithgroup.congroup_gamma1 import Gamma1_class
from sage.modular.dirichlet import DirichletGroup
from sage.modules.module import Module
from sage.rings.fast_arith import prime_range
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.primes import Primes
from sage.structure.richcmp import richcmp_method, richcmp

from .finite_subgroup import FiniteSubgroup


@richcmp_method
class RationalTorsionSubgroup(FiniteSubgroup):
    """
    The torsion subgroup of a modular abelian variety.
    """
    def __init__(self, abvar):
        """
        Create the torsion subgroup.

        INPUT:

        - ``abvar`` -- a modular abelian variety

        EXAMPLES::

            sage: T = J0(14).rational_torsion_subgroup(); T
            Torsion subgroup of Abelian variety J0(14) of dimension 1
            sage: type(T)
            <class 'sage.modular.abvar.torsion_subgroup.RationalTorsionSubgroup_with_category'>
        """
        FiniteSubgroup.__init__(self, abvar)

    def _repr_(self):
        """
        Return string representation of this torsion subgroup.

        EXAMPLES::

            sage: T = J1(13).rational_torsion_subgroup(); T
            Torsion subgroup of Abelian variety J1(13) of dimension 2
            sage: T._repr_()
            'Torsion subgroup of Abelian variety J1(13) of dimension 2'
        """
        return "Torsion subgroup of %s" % self.abelian_variety()

    def __richcmp__(self, other, op):
        """
        Compare torsion subgroups.

        INPUT:

        - ``other`` -- an object

        If other is a torsion subgroup, the abelian varieties are compared.
        Otherwise, the generic behavior for finite abelian variety
        subgroups is used.

        EXAMPLES::

            sage: G = J0(11).rational_torsion_subgroup(); H = J0(13).rational_torsion_subgroup()
            sage: G == G
            True
            sage: G < H   # since 11 < 13
            True
            sage: G > H
            False
        """
        if isinstance(other, RationalTorsionSubgroup):
            return richcmp(self.abelian_variety(), other.abelian_variety(), op)
        return FiniteSubgroup.__richcmp__(self, other, op)

    def order(self, proof=True):
        """
        Return the order of the torsion subgroup of this modular abelian
        variety.

        This function may fail if the multiple obtained by counting points
        modulo `p` exceeds the divisor obtained from the rational cuspidal
        subgroup.

        The computation of the rational torsion order of J1(p) is conjectural
        and will only be used if ``proof=False``. See Section 6.2.3 of [CES2003]_.

        INPUT:

        - ``proof`` -- boolean (default: ``True``)

        OUTPUT: the order of this torsion subgroup

        EXAMPLES::

            sage: A = J0(11)
            sage: A.rational_torsion_subgroup().order()
            5
            sage: A = J0(23)
            sage: A.rational_torsion_subgroup().order()
            11
            sage: T = J0(37)[1].rational_torsion_subgroup()
            sage: T.order()
            3

            sage: J = J1(13)
            sage: J.rational_torsion_subgroup().order()
            19

        Sometimes the order can only be computed with ``proof=False``. ::

            sage: J = J1(23)
            sage: J.rational_torsion_subgroup().order()
            Traceback (most recent call last):
            ...
            RuntimeError: Unable to compute order of torsion subgroup
            (it is in [408991, 9406793])

            sage: J.rational_torsion_subgroup().order(proof=False)
            408991
        """
        O = self.possible_orders(proof=proof)
        if len(O) == 1:
            n = O[0]
            self._order = n
            return n
        raise RuntimeError("Unable to compute order of torsion subgroup (it is in %s)" % O)

    def lattice(self):
        """
        Return lattice that defines this torsion subgroup, if possible.

        .. warning::

           There is no known algorithm in general to compute the
           rational torsion subgroup. Use rational_cusp_group to
           obtain a subgroup of the rational torsion subgroup in
           general.

        EXAMPLES::

            sage: J0(11).rational_torsion_subgroup().lattice()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [  1   0]
            [  0 1/5]

        The following fails because in fact I know of no (reasonable)
        algorithm to provably compute the torsion subgroup in general.

        ::

            sage: T = J0(33).rational_torsion_subgroup()
            sage: T.lattice()
            Traceback (most recent call last):
            ...
            NotImplementedError: unable to compute the rational torsion subgroup
            in this case (there is no known general algorithm yet)

        The problem is that the multiple of the order obtained by counting
        points over finite fields is twice the divisor of the order got
        from the rational cuspidal subgroup.

        ::

            sage: T.multiple_of_order(30)
            200
            sage: J0(33).rational_cusp_subgroup().order()
            100
        """
        A = self.abelian_variety()
        if A.dimension() == 0:
            return []
        R = A.rational_cusp_subgroup()
        if R.order() == self.multiple_of_order():
            return R.lattice()
        else:
            raise NotImplementedError("unable to compute the rational torsion subgroup in this case (there is no known general algorithm yet)")

    def possible_orders(self, proof=True):
        """
        Return the possible orders of this torsion subgroup. Outside of special
        cases, this is done by computing a divisor and multiple of the order.

        INPUT:

        - ``proof`` -- boolean (default: ``True``)

        OUTPUT: an array of positive integers

        The computation of the rational torsion order of J1(p) is conjectural
        and will only be used if ``proof=False``. See Section 6.2.3 of [CES2003]_.

        EXAMPLES::

            sage: J0(11).rational_torsion_subgroup().possible_orders()
            [5]
            sage: J0(33).rational_torsion_subgroup().possible_orders()
            [100, 200]

            sage: J1(13).rational_torsion_subgroup().possible_orders()
            [19]
            sage: J1(16).rational_torsion_subgroup().possible_orders()
            [1, 2, 4, 5, 10, 20]
        """
        try:
            if proof:
                return self._possible_orders
            else:
                return self._possible_orders_proof_false
        except AttributeError:
            pass

        A = self.abelian_variety()
        N = A.level()
        # return the order of the cuspidal subgroup in the J0(p) case
        if A.is_J0() and N.is_prime():
            self._possible_orders = [QQ((A.level()-1)/12).numerator()]
            self._possible_orders_proof_false = self._possible_orders
            return self._possible_orders

        # the elliptic curve case
        if A.dimension() == 1:
            self._possible_orders = [A.elliptic_curve().torsion_order()]
            self._possible_orders_proof_false = self._possible_orders
            return self._possible_orders

        # the conjectural J1(p) case
        if not proof and A.is_J1() and N.is_prime():
            epsilons = [epsilon for epsilon in DirichletGroup(N)
                        if not epsilon.is_trivial() and epsilon.is_even()]
            bernoullis = [epsilon.bernoulli(2) for epsilon in epsilons]
            self._possible_orders_proof_false = [ZZ(N/(2**(N-3))*prod(bernoullis))]
            return self._possible_orders_proof_false

        u = self.multiple_of_order()
        l = self.divisor_of_order()

        assert u % l == 0
        O = [l * d for d in divisors(u//l)]
        self._possible_orders = O
        if u == l:
            self._possible_orders_proof_false = O
        return O

    def divisor_of_order(self):
        """
        Return a divisor of the order of this torsion subgroup of a modular
        abelian variety.

        OUTPUT: a divisor of this torsion subgroup

        EXAMPLES::

            sage: t = J0(37)[1].rational_torsion_subgroup()
            sage: t.divisor_of_order()
            3

            sage: J = J1(19)
            sage: J.rational_torsion_subgroup().divisor_of_order()
            4383

            sage: J = J0(45)
            sage: J.rational_cusp_subgroup().order()
            32
            sage: J.rational_cuspidal_subgroup().order()
            64
            sage: J.rational_torsion_subgroup().divisor_of_order()
            64
        """
        try:
            return self._divisor_of_order
        except AttributeError:
            pass

        A = self.abelian_variety()
        N = A.level()

        if A.dimension() == 0:
            self._divisor_of_order = ZZ(1)
            return self._divisor_of_order

        # return the order of the cuspidal subgroup in the J0(p) case
        if A.is_J0() and N.is_prime():
            self._divisor_of_order = QQ((A.level()-1)/12).numerator()
            return self._divisor_of_order

        # The elliptic curve case
        if A.dimension() == 1:
            self._divisor_of_order = A.elliptic_curve().torsion_order()
            return self._divisor_of_order

        # The J1(p) case
        if A.is_J1() and N.is_prime():
            epsilons = [epsilon for epsilon in DirichletGroup(N)
                        if not epsilon.is_trivial() and epsilon.is_even()]
            bernoullis = [epsilon.bernoulli(2) for epsilon in epsilons]
            self._divisor_of_order = ZZ(N/(2**(N-3))*prod(bernoullis))
            return self._divisor_of_order

        # The Gamma0 case
        if all(isinstance(G, Gamma0_class) for G in A.groups()):
            self._divisor_of_order = A.rational_cuspidal_subgroup().order()
            return self._divisor_of_order

        # Unhandled case
        self._divisor_of_order = ZZ(1)
        return self._divisor_of_order

    def multiple_of_order(self, maxp=None, proof=True):
        """
        Return a multiple of the order.

        INPUT:

        - ``proof`` -- boolean (default: ``True``)

        The computation of the rational torsion order of J1(p) is conjectural
        and will only be used if proof=False. See Section 6.2.3 of [CES2003]_.

        EXAMPLES::

            sage: J = J1(11); J
            Abelian variety J1(11) of dimension 1
            sage: J.rational_torsion_subgroup().multiple_of_order()
            5

            sage: J = J0(17)
            sage: J.rational_torsion_subgroup().order()
            4

        This is an example where proof=False leads to a better bound and better
        performance. ::

            sage: J = J1(23)
            sage: J.rational_torsion_subgroup().multiple_of_order()  # long time (2s)
            9406793
            sage: J.rational_torsion_subgroup().multiple_of_order(proof=False)
            408991
        """

        try:
            if proof:
                return self._multiple_of_order
            else:
                return self._multiple_of_order_proof_false
        except AttributeError:
            pass

        A = self.abelian_variety()
        N = A.level()

        if A.dimension() == 0:
            self._multiple_of_order = ZZ(1)
            self._multiple_of_order_proof_false = self._multiple_of_order
            return self._multiple_of_order

        # return the order of the cuspidal subgroup in the J0(p) case
        if A.is_J0() and N.is_prime():
            self._multiple_of_order = QQ((A.level()-1)/12).numerator()
            self._multiple_of_order_proof_false = self._multiple_of_order
            return self._multiple_of_order

        # The elliptic curve case
        if A.dimension() == 1:
            self._multiple_of_order = A.elliptic_curve().torsion_order()
            self._multiple_of_order_proof_false = self._multiple_of_order
            return self._multiple_of_order

        # The conjectural J1(p) case
        if not proof and A.is_J1() and N.is_prime():
            epsilons = [epsilon for epsilon in DirichletGroup(N)
                        if not epsilon.is_trivial() and epsilon.is_even()]
            bernoullis = [epsilon.bernoulli(2) for epsilon in epsilons]
            self._multiple_of_order_proof_false = ZZ(N/(2**(N-3))*prod(bernoullis))
            return self._multiple_of_order_proof_false

        # The Gamma0 and Gamma1 case
        if all(isinstance(G, (Gamma0_class, Gamma1_class)) for G in A.groups()):
            self._multiple_of_order = self.multiple_of_order_using_frobp()
            return self._multiple_of_order

        # Unhandled case
        raise NotImplementedError("No implemented algorithm")

    def multiple_of_order_using_frobp(self, maxp=None):
        """
        Return a multiple of the order of this torsion group.

        In the `Gamma_0` case, the multiple is computed using characteristic
        polynomials of Hecke operators of odd index not dividing the level. In
        the `Gamma_1` case, the multiple is computed by expressing the
        frobenius polynomial in terms of the characteristic polynomial of left
        multiplication by `a_p` for odd primes p not dividing the level.

        INPUT:

        - ``maxp`` -- (default: ``None``) if ``maxp`` is ``None``, return gcd
          of best bound computed so far with bound obtained by computing GCD's
          of orders modulo `p` until this gcd stabilizes for 3 successive
          primes. If ``maxp`` is given, just use all primes up to and including
          ``maxp``.

        EXAMPLES::

            sage: J = J0(11)
            sage: G = J.rational_torsion_subgroup()
            sage: G.multiple_of_order_using_frobp(11)
            5

        Increasing maxp may yield a tighter bound. If maxp=None, then Sage
        will use more primes until the multiple stabilizes for 3 successive
        primes.  ::

            sage: J = J0(389)
            sage: G = J.rational_torsion_subgroup(); G
            Torsion subgroup of Abelian variety J0(389) of dimension 32
            sage: G.multiple_of_order_using_frobp()
            97
            sage: [G.multiple_of_order_using_frobp(p) for p in prime_range(3,11)]
            [92645296242160800, 7275, 291]
            sage: [G.multiple_of_order_using_frobp(p) for p in prime_range(3,13)]
            [92645296242160800, 7275, 291, 97]
            sage: [G.multiple_of_order_using_frobp(p) for p in prime_range(3,19)]
            [92645296242160800, 7275, 291, 97, 97, 97]

        We can compute the multiple of order of the torsion subgroup for Gamma0
        and Gamma1 varieties, and their products. ::

            sage: A = J0(11) * J0(33)
            sage: A.rational_torsion_subgroup().multiple_of_order_using_frobp()
            1000

            sage: A = J1(23)
            sage: A.rational_torsion_subgroup().multiple_of_order_using_frobp()
            9406793
            sage: A.rational_torsion_subgroup().multiple_of_order_using_frobp(maxp=50)
            408991

            sage: A = J1(19) * J0(21)
            sage: A.rational_torsion_subgroup().multiple_of_order_using_frobp()
            35064

        The next example illustrates calling this function with a larger
        input and how the result may be cached when maxp is None::

            sage: T = J0(43)[1].rational_torsion_subgroup()
            sage: T.multiple_of_order_using_frobp()
            14
            sage: T.multiple_of_order_using_frobp(50)
            7
            sage: T.multiple_of_order_using_frobp()
            7

        This function is not implemented for general congruence subgroups
        unless the dimension is zero. ::

            sage: A = JH(13,[2]); A
            Abelian variety J0(13) of dimension 0
            sage: A.rational_torsion_subgroup().multiple_of_order_using_frobp()
            1

            sage: A = JH(15, [2]); A
            Abelian variety JH(15,[2]) of dimension 1
            sage: A.rational_torsion_subgroup().multiple_of_order_using_frobp()
            Traceback (most recent call last):
            ...
            NotImplementedError: torsion multiple only implemented for Gamma0 and Gamma1
        """
        if maxp is None:
            try:
                return self.__multiple_of_order_using_frobp
            except AttributeError:
                pass
        A = self.abelian_variety()
        if A.dimension() == 0:
            T = ZZ.one()
            self.__multiple_of_order_using_frobp = T
            return T
        if not all(isinstance(G, (Gamma0_class, Gamma1_class)) for G in A.groups()):
            raise NotImplementedError("torsion multiple only implemented for Gamma0 and Gamma1")

        bnd = ZZ.zero()
        N = A.level()
        cnt = 0
        if maxp is None:
            X = Primes()
        else:
            X = prime_range(maxp+1)
        for p in X:
            if (2*N) % p == 0:
                continue

            if (len(A.groups()) == 1 and isinstance(A.groups()[0], Gamma0_class)):
                f = A.hecke_polynomial(p)
                b = ZZ(f(p+1))
            else:
                from .constructor import AbelianVariety
                D = [AbelianVariety(f) for f in
                     A.newform_decomposition('a')]
                b = 1
                for simple in D:
                    G = simple.newform_level()[1]
                    if isinstance(G, Gamma0_class):
                        f = simple.hecke_polynomial(p)
                        b *= ZZ(f(p+1))
                    else:
                        f = simple.newform('a')
                        Kf = f.base_ring()
                        eps = f.character()
                        Qe = eps.base_ring()

                        if Kf != QQ:
                            # relativize number fields to compute charpoly of
                            # left multiplication of ap on Kf as a Qe-vector
                            # space.
                            Lf = Kf.relativize(Qe.gen(), 'a')
                            to_Lf = Lf.structure()[1]

                            name = Kf._names[0]
                            ap = to_Lf(f.modular_symbols(1).eigenvalue(p, name))

                            G_ps = ap.matrix().charpoly()
                            b *= ZZ(Qe(G_ps(1 + to_Lf(eps(p))*p)).norm())
                        else:
                            ap = f.modular_symbols(1).eigenvalue(p)
                            b *= ZZ(1 + eps(p)*p - ap)

            if bnd == 0:
                bnd = b
            else:
                bnd_last = bnd
                bnd = ZZ(gcd(bnd, b))
                if bnd == bnd_last:
                    cnt += 1
                else:
                    cnt = 0
                if maxp is None and cnt >= 2:
                    break

        # The code below caches the computed bound and
        # will be used if this function is called
        # again with maxp equal to None (the default).
        if maxp is None:
            # maxp is None but self.__multiple_of_order_using_frobp  is
            # not set, since otherwise we would have immediately
            # returned at the top of this function
            self.__multiple_of_order_using_frobp = bnd
        else:
            # maxp is given -- record new info we get as
            # a gcd...
            try:
                self.__multiple_of_order_using_frobp = \
                        gcd(self.__multiple_of_order_using_frobp, bnd)
            except AttributeError:
                # ... except in the case when
                # self.__multiple_of_order_using_frobp was never set.  In this
                # case, we just set it as long as the gcd stabilized for 3 in a
                # row.
                if cnt >= 2:
                    self.__multiple_of_order_using_frobp = bnd
        return bnd


class QQbarTorsionSubgroup(Module):

    Element = TorsionPoint

    def __init__(self, abvar):
        """
        Group of all torsion points over the algebraic closure on an
        abelian variety.

        INPUT:

        - ``abvar`` -- an abelian variety

        EXAMPLES::

            sage: A = J0(23)
            sage: A.qbar_torsion_subgroup()                                             # needs sage.rings.number_field
            Group of all torsion points in QQbar on Abelian variety J0(23) of dimension 2
        """
        self.__abvar = abvar
        Module.__init__(self, ZZ)

    def _repr_(self):
        """
        Print representation of QQbar points.

        OUTPUT: string

        EXAMPLES::

            sage: J0(23).qbar_torsion_subgroup()._repr_()                               # needs sage.rings.number_field
            'Group of all torsion points in QQbar on Abelian variety J0(23) of dimension 2'
        """
        return 'Group of all torsion points in QQbar on %s' % self.__abvar

    def field_of_definition(self):
        """
        Return the field of definition of this subgroup. Since this is the
        group of all torsion it is defined over the base field of this
        abelian variety.

        OUTPUT: a field

        EXAMPLES::

            sage: J0(23).qbar_torsion_subgroup().field_of_definition()                  # needs sage.rings.number_field
            Rational Field
        """
        return self.__abvar.base_field()

    def _element_constructor_(self, x):
        r"""
        Create an element in this torsion subgroup.

        INPUT:

        - ``x`` -- vector in `\QQ^{2d}`

        OUTPUT: torsion point

        EXAMPLES::

            sage: P = J0(23).qbar_torsion_subgroup()([1,1/2,3/4,2]); P                  # needs sage.rings.number_field
            [(1, 1/2, 3/4, 2)]
            sage: P.order()                                                             # needs sage.rings.number_field
            4
        """
        v = self.__abvar.vector_space()(x)
        return self.element_class(self, v)

    def abelian_variety(self):
        """
        Return the abelian variety that this is the set of all torsion
        points on.

        OUTPUT: abelian variety

        EXAMPLES::

            sage: J0(23).qbar_torsion_subgroup().abelian_variety()                      # needs sage.rings.number_field
            Abelian variety J0(23) of dimension 2
        """
        return self.__abvar
