# sage_setup: distribution = sagemath-polyhedra
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
