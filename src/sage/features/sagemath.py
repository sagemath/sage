r"""
Check for SageMath Python modules
"""
from . import PythonModule


class sage__combinat(PythonModule):

    def __init__(self):
        # sage.combinat will be a namespace package.
        # Testing whether sage.combinat itself can be imported is meaningless.
        # Hence, we test a Python module within the package.
        PythonModule.__init__(self, 'sage.combinat.combinations')


class sage__graphs(PythonModule):

    def __init__(self):
        PythonModule.__init__(self, 'sage.graphs.graph')


class sage__rings__number_field(PythonModule):

    def __init__(self):
        PythonModule.__init__(self, 'sage.rings.number_field_element')


class sage__rings__real_double(PythonModule):

    def __init__(self):
        PythonModule.__init__(self, 'sage.rings.real_double')


class sage__symbolic(PythonModule):

    def __init__(self):
        PythonModule.__init__(self, 'sage.symbolic.expression', spkg="sagemath_symbolics")


def sage_optional_tags():
    """
    Return tags for conditionalizing doctests.

    These tags are named after Python packages/modules (e.g., :mod:`~sage.symbolic`),
    not distribution packages (``sagemath-symbolics``).

    This is motivated by a separation of concerns: The author of a module that depends
    on some functionality provided by a Python module usually already knows the
    name of the Python module, so we do not want to force the author to also
    know about the distribution package that provides the Python module.

    Instead, we associate distribution packages to Python modules in
    :mod:`sage.features.sagemath` via the ``spkg`` parameter of :class:`PythonModule``.
    """
    if sage__combinat().is_present():
        yield 'sage.combinat'
    if sage__graphs().is_present():
        yield 'sage.graphs'
    if sage__rings__number_field().is_present():
        yield 'sage.rings.number_field'
    if sage__rings__real_double().is_present():
        yield 'sage.rings.real_double'
    if sage__symbolic().is_present():
        yield 'sage.symbolic'
