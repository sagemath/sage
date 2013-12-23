"""
The global ``dev`` object in Sage

EXAMPLES::

    sage: dev
    SageDev()
    sage: type(dev._get_object())
    <class 'sage.dev.test.sagedev.DoctestSageDevWrapper'>
"""


from sage.doctest import DOCTEST_MODE
if DOCTEST_MODE:
    from sage.dev.test.sagedev import DoctestSageDevWrapper
    from sage.dev.test.config import DoctestConfig
    from sage.dev.test.trac_server import DoctestTracServer
    dev = DoctestSageDevWrapper(DoctestConfig(), DoctestTracServer())
else:
    from sagedev import SageDev
    from sagedev_wrapper import SageDevWrapper
    dev = SageDevWrapper(SageDev())
