from sage.misc.lazy_import import lazy_import
lazy_import('sage.manifolds.manifold', 'Manifold')
lazy_import('sage.manifolds.differentiable.examples.real_line', ('OpenInterval', 'RealLine'),
            deprecation=31881)
lazy_import('sage.manifolds.differentiable.examples.euclidean', 'EuclideanSpace')
lazy_import('sage.manifolds', 'catalog', 'manifolds')
