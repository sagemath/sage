import pytest
from sage.numerical.backends.generic_backend import GenericBackend
from sage.structure.sage_object import SageObject
from sage.structure.sage_object_test import SageObjectTests


class GenericBackendTests(SageObjectTests):

    @pytest.fixture
    def backend(self, *args, **kwargs) -> GenericBackend:
        raise NotImplementedError

    @pytest.fixture
    def sage_object(self, backend) -> SageObject:
        return backend

    def test_ncols_nonnegative(self, backend: GenericBackend):
        assert backend.ncols() >= 0

    def test_sage_unittest_testsuite(self, sage_object: SageObject, run_test_suite):
        # TODO: Remove this test as soon as all old test methods are migrated
        run_test_suite(sage_object, verbose=True, skip="_test_pickling")
