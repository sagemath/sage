from __future__ import annotations

import pytest

from sage.numerical.backends.generic_backend_test import GenericBackendTests

CoinBackend = pytest.importorskip('sage.numerical.backends.coin_backend').CoinBackend


class TestCoinBackend(GenericBackendTests):

    @pytest.fixture
    def backend(self):
        return CoinBackend()

    def test_sage_unittest_testsuite(self, sage_object):
        from sage.misc.sage_unittest import TestSuite

        # Keep the exclusions used by CoinBackend's existing doctest.
        TestSuite(sage_object).run(raise_on_failure=True,
                                  skip=['_test_pickling', '_test_solve'])

    def test_solver_selection(self):
        from sage.numerical.mip import MixedIntegerLinearProgram

        program = MixedIntegerLinearProgram(solver='Coin')
        assert isinstance(program.get_backend(), CoinBackend)

    def test_feature(self):
        from sage.features.mip_backends import COIN

        assert COIN().name == 'cbc'
        assert COIN().spkg == 'cbc'
        assert COIN().is_present()
