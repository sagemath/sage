# sage_setup: distribution = sagemath-objects

import pytest
from sage.structure.sage_object import SageObject


class SageObjectTests:

    @pytest.fixture
    def sage_object(self, *args, **kwargs) -> SageObject:
        raise NotImplementedError

    def test_sage_unittest_testsuite(self, sage_object: SageObject):
        """
        Subclasses should override this method if they need to skip some tests.
        """
        # TODO: Remove this test as soon as all old test methods are migrated
        from sage.misc.sage_unittest import TestSuite
        TestSuite(sage_object).run(verbose=True, raise_on_failure=True)
