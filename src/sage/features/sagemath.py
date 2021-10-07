r"""
Check for SageMath Python modules
"""
from . import PythonModule
from .join_feature import JoinFeature


class sage__combinat(JoinFeature):

    def __init__(self):
        # sage.combinat will be a namespace package.
        # Testing whether sage.combinat itself can be imported is meaningless.
        # Hence, we test a Python module within the package.
        JoinFeature.__init__(self, 'sage.combinat',
                             [PythonModule('sage.combinat.combinations')])


class sage__graphs(JoinFeature):

    def __init__(self):
        JoinFeature.__init__(self, 'sage.graphs',
                             [PythonModule('sage.graphs.graph')])


class sage__graphs__bliss(PythonModule):

    def __init__(self):
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_bliss' later
        PythonModule.__init__(self, 'sage.graphs.bliss', spkg='bliss')


class sage__graphs__graph_decompositions__tdlib(PythonModule):

    def __init__(self):
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_tdlib' later
        PythonModule.__init__(self, 'sage.graphs.graph_decompositions.tdlib', spkg='tdlib')


class sage__graphs__mcqd(PythonModule):

    def __init__(self):
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_mcqd' later
        PythonModule.__init__(self, 'sage.graphs.mcqd', spkg='mcqd')


class sage__matrix__matrix_gfpn_dense(PythonModule):

    def __init__(self):
        # Currently part of sagemath_standard, conditionally built.
        # Will be changed to spkg='sagemath_meataxe' later
        PythonModule.__init__(self, 'sage.matrix.matrix_gfpn_dense', spkg='meataxe')


class sage__plot(JoinFeature):

    def __init__(self):
        JoinFeature.__init__(self, 'sage.plot',
                             [PythonModule('sage.plot.plot')])


class sage__rings__number_field(JoinFeature):

    def __init__(self):
        JoinFeature.__init__(self, 'sage.rings.number_field',
                             [PythonModule('sage.rings.number_field.number_field_element')])


class sage__rings__real_double(PythonModule):

    def __init__(self):
        PythonModule.__init__(self, 'sage.rings.real_double')


class sage__symbolic(JoinFeature):

    def __init__(self):
        JoinFeature.__init__(self, 'sage.symbolic',
                             [PythonModule('sage.symbolic.expression')],
                             spkg="sagemath_symbolics")


def sage_features():
    """
    Return tags for conditionalizing doctests.

    These tags are named after Python packages/modules (e.g., :mod:`~sage.symbolic`),
    not distribution packages (``sagemath-symbolics``).

    This design is motivated by a separation of concerns: The author of a module that depends
    on some functionality provided by a Python module usually already knows the
    name of the Python module, so we do not want to force the author to also
    know about the distribution package that provides the Python module.

    Instead, we associate distribution packages to Python modules in
    :mod:`sage.features.sagemath` via the ``spkg`` parameter of :class:`Feature`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage_features
        sage: list(sage_features())  # random
        [Feature('sage.graphs'),
         Feature('sage.graphs.bliss'),
         Feature('sage.plot'),
         Feature('sage.rings.real_double')]
    """
    for feature in [sage__combinat(),
                    sage__graphs(),
                    sage__graphs__bliss(),
                    sage__graphs__graph_decompositions__tdlib(),
                    sage__graphs__mcqd(),
                    sage__matrix__matrix_gfpn_dense(),
                    sage__plot(),
                    sage__rings__number_field(),
                    sage__rings__real_double(),
                    sage__symbolic()]:
        if feature.is_present():
            yield feature
