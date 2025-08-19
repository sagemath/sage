# sage.doctest: needs sage.groups
r"""
Testing arithmetic subgroup
"""
################################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
#
################################################################################
from __future__ import annotations

from .arithgroup_perm import ArithmeticSubgroup_Permutation, EvenArithmeticSubgroup_Permutation, OddArithmeticSubgroup_Permutation
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor as Gamma0
from sage.modular.arithgroup.congroup_gamma1 import Gamma1_constructor as Gamma1
from sage.modular.arithgroup.congroup_gammaH import GammaH_constructor as GammaH
from sage.modular.arithgroup.congroup_gamma import Gamma_constructor as Gamma
from sage.rings.finite_rings.integer_mod_ring import Zmod

import sage.misc.prandom as prandom
from sage.misc.timing import cputime


def random_even_arithgroup(index, nu2_max=None, nu3_max=None):
    r"""
    Return a random even arithmetic subgroup.

    EXAMPLES::

        sage: import sage.modular.arithgroup.tests as tests
        sage: G = tests.random_even_arithgroup(30); G # random
        Arithmetic subgroup of index 30
        sage: G.is_even()
        True
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup

    test = False

    if nu2_max is None:
        nu2_max = index // 5
    elif nu2_max == 0:
        assert index % 2 == 0
    if nu3_max is None:
        nu3_max = index // 7
    elif nu3_max == 0:
        assert index % 3 == 0

    while not test:
        nu2 = prandom.randint(0, nu2_max)
        nu2 = index % 2 + nu2 * 2
        nu3 = prandom.randint(0, nu3_max)
        nu3 = index % 3 + nu3 * 3

        l = list(range(1, index + 1))
        prandom.shuffle(l)
        S2: list[tuple] = []
        for i in range(nu2):
            S2.append((l[i],))
        for i in range(nu2, index, 2):
            S2.append((l[i], l[i + 1]))
        prandom.shuffle(l)
        S3: list[tuple] = []
        for i in range(nu3):
            S3.append((l[i],))
        for i in range(nu3, index, 3):
            S3.append((l[i], l[i + 1], l[i + 2]))
        G = PermutationGroup([S2, S3])
        test = G.is_transitive()

    return ArithmeticSubgroup_Permutation(S2=S2, S3=S3)


def random_odd_arithgroup(index, nu3_max=None):
    r"""
    Return a random odd arithmetic subgroup.

    EXAMPLES::

        sage: from sage.modular.arithgroup.tests import random_odd_arithgroup
        sage: G = random_odd_arithgroup(20); G #random
        Arithmetic subgroup of index 20
        sage: G.is_odd()
        True
    """
    assert index % 4 == 0
    G = random_even_arithgroup(index // 2, nu2_max=0, nu3_max=nu3_max)
    return G.one_odd_subgroup(random=True)


class Test:
    r"""
    Testing class for arithmetic subgroup implemented via permutations.
    """

    def __init__(self, index=20, index_max=50, odd_probability=0.5):
        r"""
        Create an arithmetic subgroup testing object.

        INPUT:

        - ``index`` -- the index of random subgroup to test

        - ``index_max`` -- the maximum index for congruence subgroup to test

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test()
            Arithmetic subgroup testing class
        """
        self.congroups = []
        i = 1
        self.odd_probability = odd_probability
        if index % 4:
            self.odd_probability = 0
        while Gamma(i).index() < index_max:
            self.congroups.append(Gamma(i))
            i += 1
        i = 1
        while Gamma0(i).index() < index_max:
            self.congroups.append(Gamma0(i))
            i += 1
        i = 2
        while Gamma1(i).index() < index_max:
            self.congroups.append(Gamma1(i))
            M = Zmod(i)
            U = [x for x in M if x.is_unit()]
            for j in range(1, len(U) - 1):
                self.congroups.append(GammaH(i, prandom.sample(U, j)))
            i += 1

        self.index = index

    def __repr__(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().__repr__()
            'Arithmetic subgroup testing class'
        """
        return "Arithmetic subgroup testing class"

    def _do(self, name):
        """
        Perform the test 'test_name', where name is specified as an
        argument. This function exists to avoid a call to eval.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test()._do("random")
            test_random
            ...
        """
        print("test_%s" % name)
        Test.__dict__["test_%s" % name](self)

    def random(self, seconds=0):
        """
        Perform random tests for a given number of seconds, or
        indefinitely if seconds is not specified.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().random(1)
            test_random
            ...
        """
        self.test("random", seconds)

    def test(self, name, seconds=0):
        """
        Repeatedly run 'test_name', where name is passed as an
        argument. If seconds is nonzero, run for that many seconds. If
        seconds is 0, run indefinitely.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: T = Test()
            sage: T.test('relabel',seconds=1)
            test_relabel
            ...
            sage: T.test('congruence_groups',seconds=1)
            test_congruence_groups
            ...
            sage: T.test('contains',seconds=1)
            test_contains
            ...
            sage: T.test('todd_coxeter',seconds=1)
            test_todd_coxeter
            ...
        """
        seconds = float(seconds)
        total = cputime()
        n = 1
        while seconds == 0 or cputime(total) < seconds:
            s = "** test_dimension: number %s" % n
            if seconds > 0:
                s += " (will stop after about %s seconds)" % seconds
            t = cputime()
            self._do(name)
            print("\ttime=%s\telapsed=%s" % (cputime(t), cputime(total)))
            n += 1

    def test_random(self):
        """
        Do a random test from all the possible tests.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().test_random() #random
            Doing random test
        """
        tests = [a for a in Test.__dict__
                 if a[:5] == "test_" and a != "test_random"]
        name = prandom.choice(tests)
        print("Doing random test %s" % name)
        Test.__dict__[name](self)

    def test_relabel(self):
        r"""
        Try the function canonical labels for a random even modular subgroup.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().test_relabel() # random
        """
        if prandom.uniform(0, 1) < self.odd_probability:
            G = random_odd_arithgroup(self.index)
        else:
            G = random_even_arithgroup(self.index)

        G.relabel()
        s2 = G._S2
        s3 = G._S3
        l = G._L
        r = G._R

        # 0 should be stabilized by the mapping
        # used for renumbering so we start at 1
        p = list(range(1, self.index))

        for _ in range(10):
            prandom.shuffle(p)
            # we add 0 to the mapping
            pp = [0] + p
            ss2 = [None] * self.index
            ss3 = [None] * self.index
            ll = [None] * self.index
            rr = [None] * self.index
            for i in range(self.index):
                ss2[pp[i]] = pp[s2[i]]
                ss3[pp[i]] = pp[s3[i]]
                ll[pp[i]] = pp[l[i]]
                rr[pp[i]] = pp[r[i]]
            if G.is_even():
                GG = EvenArithmeticSubgroup_Permutation(ss2, ss3, ll, rr)
            else:
                GG = OddArithmeticSubgroup_Permutation(ss2, ss3, ll, rr)
            GG.relabel()

            for elt in ['_S2', '_S3', '_L', '_R']:
                if getattr(G, elt) != getattr(GG, elt):
                    print("s2 = %s" % str(s2))
                    print("s3 = %s" % str(s3))
                    print("ss2 = %s" % str(ss2))
                    print("ss3 = %s" % str(ss3))
                    print("pp = %s" % str(pp))
                    raise AssertionError("%s does not coincide" % elt)

    def test_congruence_groups(self):
        r"""
        Check whether the different implementations of methods for congruence
        groups and generic arithmetic group by permutations return the same
        results.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().test_congruence_groups() #random
        """
        G = prandom.choice(self.congroups)
        GG = G.as_permutation_group()

        if not GG.is_congruence():
            raise AssertionError("Hsu congruence test failed")

        methods = [
            'index',
            'is_odd',
            'is_even',
            'is_normal',
            'ncusps',
            'nregcusps',
            'nirregcusps',
            'nu2',
            'nu3',
            'generalised_level']

        for f in methods:
            if getattr(G, f)() != getattr(GG, f)():
                raise AssertionError("results of %s does not coincide for %s" % (f, G))

        if sorted(G.cusp_width(c) for c in G.cusps()) != GG.cusp_widths():
            raise AssertionError("Cusps widths are different for %s" % G)

        for _ in range(20):
            m = GG.random_element()
            if m not in G:
                raise AssertionError("random element generated by perm. group not in %s" % str(G))

    def test_contains(self):
        r"""
        Test whether the random generator for arithgroup perms gives matrices in
        the group.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().test_contains() #random
        """
        if prandom.uniform(0, 1) < self.odd_probability:
            G = random_odd_arithgroup(self.index)
        else:
            G = random_even_arithgroup(self.index)

        for _ in range(20):
            g = G.random_element()
            if G.random_element() not in G:
                raise AssertionError("%s not in %s" % (g, G))

    def test_spanning_trees(self):
        r"""
        Test coset representatives obtained from spanning trees for even
        subgroup (Kulkarni's method with generators ``S2``, ``S3`` and Verrill's
        method with generators ``L``, ``S2``).

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().test_spanning_trees() #random
        """
        from sage.misc.misc_c import prod
        from .all import SL2Z
        from .arithgroup_perm import S2m, S3m, Lm

        G = random_even_arithgroup(self.index)

        m = {'l': Lm, 's': S2m}
        tree, reps, wreps, gens = G._spanning_tree_verrill()
        assert reps[0] == SL2Z([1, 0, 0, 1])
        assert wreps[0] == ''
        for i in range(1, self.index):
            assert prod(m[letter] for letter in wreps[i]) == reps[i]
        tree, reps, wreps, gens = G._spanning_tree_verrill(on_right=False)
        assert reps[0] == SL2Z([1, 0, 0, 1])
        assert wreps[0] == ''
        for i in range(1, self.index):
            assert prod(m[letter] for letter in wreps[i]) == reps[i]

        m = {'s2': S2m, 's3': S3m}
        tree, reps, wreps, gens = G._spanning_tree_kulkarni()
        assert reps[0] == SL2Z([1, 0, 0, 1])
        assert wreps[0] == []
        for i in range(1, self.index):
            assert prod(m[letter] for letter in wreps[i]) == reps[i]
        tree, reps, wreps, gens = G._spanning_tree_kulkarni(on_right=False)
        assert reps[0] == SL2Z([1, 0, 0, 1])
        assert wreps[0] == []
        for i in range(1, self.index):
            assert prod(m[letter] for letter in wreps[i]) == reps[i]

    def test_todd_coxeter(self):
        r"""
        Test representatives of Todd-Coxeter algorithm.

        EXAMPLES::

            sage: from sage.modular.arithgroup.tests import Test
            sage: Test().test_todd_coxeter() #random
        """
        from .all import SL2Z
        from .arithgroup_perm import S2m, S3m, Lm

        G = random_even_arithgroup(self.index)

        reps, gens, l, s2 = G.todd_coxeter_l_s2()
        assert reps[0] == SL2Z([1, 0, 0, 1])
        assert len(reps) == G.index()
        for i in range(1, len(reps)):
            assert reps[i] not in G
            assert reps[i] * S2m * ~reps[s2[i]] in G
            assert reps[i] * Lm * ~reps[l[i]] in G
            for j in range(i + 1, len(reps)):
                assert reps[i] * ~reps[j] not in G
                assert reps[j] * ~reps[i] not in G

        reps, gens, s2, s3 = G.todd_coxeter_s2_s3()
        assert reps[0] == SL2Z([1, 0, 0, 1])
        assert len(reps) == G.index()
        for i in range(1, len(reps)):
            assert reps[i] not in G
            assert reps[i] * S2m * ~reps[s2[i]] in G
            assert reps[i] * S3m * ~reps[s3[i]] in G
            for j in range(i + 1, len(reps)):
                assert reps[i] * ~reps[j] not in G
                assert reps[j] * ~reps[i] not in G
