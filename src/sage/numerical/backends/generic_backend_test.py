from __future__ import annotations
import pytest
from sage.structure.sage_object_test import SageObjectTests
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sage.structure.sage_object import SageObject
    from sage.numerical.backends.generic_backend import GenericBackend


class GenericBackendTests(SageObjectTests):

    @pytest.fixture
    def backend(self, *args, **kwargs) -> GenericBackend:
        raise NotImplementedError

    @pytest.fixture
    def sage_object(self, backend) -> SageObject:
        return backend

    def test_ncols_nonnegative(self, backend: GenericBackend):
        assert backend.ncols() >= 0

    def test_sage_unittest_testsuite(self, sage_object: SageObject):
        # TODO: Remove this test as soon as all old test methods are migrated
        from sage.misc.sage_unittest import TestSuite
        TestSuite(sage_object).run(verbose=True, raise_on_failure=True, skip='_test_pickling')
