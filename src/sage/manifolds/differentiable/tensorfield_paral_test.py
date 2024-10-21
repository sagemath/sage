# sage_setup: distribution = sagemath-symbolics
# pylint: disable=missing-function-docstring,missing-class-docstring
import pytest

from sage.manifolds.differentiable.examples.euclidean import EuclideanSpace
from sage.manifolds.differentiable.manifold import DifferentiableManifold


class TestR3VectorSpace:
    @pytest.fixture
    def manifold(self):
        return EuclideanSpace(3)

    def test_trace_using_metric_works(self, manifold: DifferentiableManifold):
        metric = manifold.metric('g')
        assert metric.trace(using=metric) == manifold.scalar_field(3)
