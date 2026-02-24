# sage.doctest: needs sage.libs.flint sage.libs.pari
"""
Testing modular symbols spaces

TESTS::

    sage: m = ModularSymbols(389)
    sage: [(g.degree(), e) for g, e in m.T(2).fcp()]
    [(1, 1), (1, 2), (2, 2), (3, 2), (6, 2), (20, 2)]
"""

# ****************************************************************************
#       Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
# ****************************************************************************

import random

from sage.misc.timing import cputime
from sage.modular import dirichlet
from sage.modular.arithgroup.congroup_gamma0 import Gamma0_constructor
from sage.modular.arithgroup.congroup_gamma1 import Gamma1_constructor
from sage.modular.modsym import modsym


class Test:
    """
    Modular symbol testing class.
    """
    def __init__(self, levels=20, weights=4, onlyg0=False, onlyg1=False,
                 onlychar=False):
        """
        Create a modular symbol testing object.

        INPUT:

        - ``levels`` -- list or integer
        - ``weights`` -- list or integer
        - ``onlyg0`` -- boolean; if ``True`` only select Gamma0 spaces for testing
        - ``onlyg1`` -- boolean; if ``True`` only select Gamma1 spaces for testing
        - ``onlychar`` -- boolean; if ``True`` only selects spaces with character for testing

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test()
            Modular symbols testing class
            sage: T = Test(weights=[3,5,7])
            sage: T.weights
            [3, 5, 7]
            sage: T = Test(levels=5) ; T.levels
            [1, 2, 3, 4, 5]
        """
        if not isinstance(levels, list):
            levels = list(range(1, int(levels) + 1))
        if not isinstance(weights, list):
            weights = list(range(2, int(weights) + 1))
        self.levels = levels
        self.weights = weights
        if not levels:
            raise RuntimeError("levels must have positive length")
        if not weights:
            raise RuntimeError("weights must have positive length")
        self.current_space = None
        self.onlyg0 = onlyg0
        self.onlyg1 = onlyg1
        self.onlychar = onlychar

    def __repr__(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test().__repr__()
            'Modular symbols testing class'
        """
        return "Modular symbols testing class"

    def _modular_symbols_space(self):
        """
        Generate a random space of modular symbols subject to
        the conditions of ``self``.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: T = Test(levels=[5],weights=[2], onlychar=True)

        Note that the sign of the generated space is always arbitrary.
            sage: T._modular_symbols_space()
            character
            level = 5, weight = 2, sign = ...
            ...
        """
        if self.onlyg0:
            which = 0
        elif self.onlyg1:
            which = 1
        elif self.onlychar:
            which = 2
        else:
            which = random.randrange(0, 3)
        if which == 0:
            print("gamma0")
            M = self._modular_symbols_space_gamma0()
        elif which == 1:
            print("gamma1")
            M = self._modular_symbols_space_gamma1()
        else:
            print("character")
            M = self._modular_symbols_space_character()
        print("\t", M)
        return M

    def _level_weight_sign(self):
        """
        Return a triple containing a random choice of level from
        self.levels, weights from self.weights, and sign chosen
        randomly from [1, 0, -1].

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test()._level_weight_sign() # random
            level = 4, weight = 3, sign = 1
            (4, 3, 1)
        """
        level = random.choice(self.levels)
        weight = random.choice(self.weights)
        sign = random.choice([-1, 0, 1])
        print("level = %s, weight = %s, sign = %s" % (level, weight, sign))
        return level, weight, sign

    def _modular_symbols_space_gamma0(self):
        """
        Return a space of modular symbols for Gamma0, with level a
        random choice from self.levels, weight from self.weights, and
        sign chosen randomly from [1, 0, -1].

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test()._modular_symbols_space_gamma0() # random
            level = 1, weight = 3, sign = 0
            Modular Symbols space of dimension 0 for Gamma_0(1) of weight 3 with sign 0 over Rational Field
        """
        level, weight, sign = self._level_weight_sign()
        M = modsym.ModularSymbols(Gamma0_constructor(level), weight, sign)
        self.current_space = M
        return M

    def _modular_symbols_space_gamma1(self):
        """
        Return a space of modular symbols for Gamma1, with level a
        random choice from self.levels, weight from self.weights, and
        sign chosen randomly from [1, 0, -1].

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test()._modular_symbols_space_gamma1() # random
            level = 3, weight = 4, sign = 0
            Modular Symbols space of dimension 2 for Gamma_1(3) of weight 4 with sign 0 over Rational Field
        """
        level, weight, sign = self._level_weight_sign()
        M = modsym.ModularSymbols(Gamma1_constructor(level), weight, sign)
        self.current_space = M
        return M

    def _modular_symbols_space_character(self):
        """
        Return a random space of modular symbols for Gamma1, with
        (possibly trivial) character, random sign, and weight a
        random choice from self.weights.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: from sage.modular.modsym.tests import Test
            sage: Test()._modular_symbols_space_character() # random
            level = 18, weight = 3, sign = 0
            Modular Symbols space of dimension 0
            and level 18, weight 3, character [1, zeta6 - 1], sign 0,
            over Cyclotomic Field of order 6 and degree 2
        """
        level, weight, sign = self._level_weight_sign()
        G = dirichlet.DirichletGroup(level)
        eps = G.random_element()
        M = modsym.ModularSymbols(eps, weight, sign)
        self.current_space = M
        return M

    def _do(self, name):
        """
        Perform the test 'test_name', where name is specified as an
        argument. This function exists to avoid a call to eval.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test()._do("random")
            test_random
            ...
        """
        print("test_%s" % name)
        Test.__dict__["test_%s" % name](self)

    #################################################################
    # The tests
    #################################################################
    def random(self, seconds=0):
        """
        Perform random tests for a given number of seconds, or
        indefinitely if seconds is not specified.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
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

            sage: from sage.modular.modsym.tests import Test
            sage: Test().test('cs_dimension',seconds=1)
            test_cs_dimension
            ...
            sage: Test().test('csnew_dimension',seconds=1)
            test_csnew_dimension
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

    def test_cs_dimension(self):
        """
        Compute the cuspidal subspace (this implicitly checks that the
        dimension is correct using formulas).

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test().test_cs_dimension() # random
            gamma0
            level = 16, weight = 3, sign = -1
            Modular Symbols space of dimension 0 for Gamma_0(16) of weight 3
            with sign -1 over Rational Field
        """
        self._modular_symbols_space().cuspidal_submodule()

    def test_csnew_dimension(self):
        """
        Compute the new cuspidal subspace and verify that the
        dimension is correct using a dimension formula.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test().test_csnew_dimension() # random
            gamma0
            level = 3, weight = 3, sign = 1
            Modular Symbols space of dimension 0 for Gamma_0(3) of weight 3
            with sign 1 over Rational Field
        """
        M = self._modular_symbols_space()
        V = M.cuspidal_submodule().new_submodule()
        d = V.dimension()
        d2 = M._cuspidal_new_submodule_dimension_formula()
        assert d == d2, \
            "Test failed for M=\"%s\", where computed dimension is %s but formula dimension is %s." % (M, d, d2)

    def test_csns_nscs(self):
        """
        Compute new cuspidal subspace in two ways and verify that the
        results are the same.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test().test_csns_nscs() # random
            gamma0
            level = 5, weight = 4, sign = 1
            Modular Symbols space of dimension 3 for Gamma_0(5) of weight 4 with sign 1 over Rational Field
        """
        M = self._modular_symbols_space()
        V1 = M.cuspidal_submodule().new_submodule()
        V2 = M.new_submodule().cuspidal_submodule()
        assert V1 == V2, "Test failed for M=\"%s\", where the new cuspidal and cuspidal new spaces are computed differently." % M
        d = M._cuspidal_new_submodule_dimension_formula()
        assert d == V1.dimension(), \
            "Test failed for M=\"%s\", where computed dimension is %s but formula dimension is %s." % (
                M, V1.dimension(), d)

    def test_decomposition(self):
        """
        Compute the decomposition of a modular symbols space, and
        verify that the sum of the dimensions of its components equals
        the dimension of the original space.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test().test_decomposition() # random
            gamma1
            level = 10, weight = 4, sign = 0
            Modular Symbols space of dimension 18 for Gamma_1(10) of weight 4 with sign 0 over Rational Field
        """
        M = self._modular_symbols_space()
        D = M.decomposition()
        assert M.dimension() == sum([A.dimension() for A in D])

    def test_dimension(self):
        """
        Compute the dimension of a modular symbols space.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test().test_dimension() # random
            gamma1
            level = 14, weight = 2, sign = -1
            Modular Symbols space of dimension 1 for Gamma_1(14) of weight 2 with sign -1 over Rational Field
        """
        self._modular_symbols_space().dimension()

    def test_random(self):
        """
        Do a random test from all the possible tests.

        EXAMPLES::

            sage: from sage.modular.modsym.tests import Test
            sage: Test().test_random() # random
            Doing random test test_csnew_dimension
            character
            level = 18, weight = 4, sign = -1
            Modular Symbols space of dimension 0 and level 18, weight 4, character [1, -1], sign -1, over Rational Field
        """
        tests = [a for a in Test.__dict__
                 if a[:5] == "test_" and a != "test_random"]
        name = random.choice(tests)
        print("Doing random test %s" % name)
        Test.__dict__[name](self)
