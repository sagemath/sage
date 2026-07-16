
import pytest
from sage.structure.sage_object import SageObject


class SageObjectTests:

    @pytest.fixture
    def sage_object(self, *args, **kwargs) -> SageObject:
        raise NotImplementedError

    def test_sage_unittest_testsuite(self, sage_object: SageObject, run_test_suite):
        """
        Subclasses should override this method if they need to skip some tests.
        """
        # TODO: Remove this test as soon as all old test methods are migrated
        run_test_suite(sage_object, verbose=True)
