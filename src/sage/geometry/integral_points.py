r"""
Cython helper methods to compute integral points in polyhedra
"""

try:
    from .integral_points_integer_dense import (
        parallelotope_points,
        ray_matrix_normal_form,
        loop_over_parallelotope_points,
        simplex_points,
        rectangular_box_points,
        print_cache,
        Inequality_generic,
        Inequality_int,
        InequalityCollection,
    )
except ImportError:
    from .integral_points_generic_dense import (
        parallelotope_points,
        ray_matrix_normal_form,
        loop_over_parallelotope_points,
        simplex_points,
        rectangular_box_points,
        print_cache,
        Inequality_generic,
        Inequality_int,
        InequalityCollection,
    )


# __all__ is needed to generate Sphinx documentation
__all__ = ['InequalityCollection', 'Inequality_generic', 'Inequality_int',
           'loop_over_parallelotope_points', 'parallelotope_points', 'print_cache',
           'ray_matrix_normal_form', 'rectangular_box_points', 'simplex_points']
