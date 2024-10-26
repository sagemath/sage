py.install_sources(
  'all.py',
  'cone.py',
  'cone_catalog.py',
  'cone_critical_angles.py',
  'convex_set.py',
  'fan.py',
  'fan_isomorphism.py',
  'fan_morphism.py',
  'hasse_diagram.py',
  'integral_points.py',
  'lattice_polytope.py',
  'linear_expression.py',
  'newton_polygon.py',
  'polyhedral_complex.py',
  'pseudolines.py',
  'relative_interior.py',
  'ribbon_graph.py',
  'toric_lattice.py',
  'toric_plotter.py',
  'voronoi_diagram.py',
  subdir: 'sage/geometry',
)

extension_data = {
  'abc' : files('abc.pyx'),
  'integral_points_generic_dense' : files('integral_points_generic_dense.pyx'),
  'integral_points_integer_dense' : files('integral_points_integer_dense.pyx'),
  'palp_normal_form' : files('palp_normal_form.pyx'),
  'point_collection' : files('point_collection.pyx'),
  'toric_lattice_element' : files('toric_lattice_element.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/geometry',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_rings],
    dependencies: [py_dep, cysignals, flint, gmp],
  )
endforeach

install_subdir('hyperbolic_space', install_dir: sage_install_dir / 'geometry')
install_subdir(
  'hyperplane_arrangement',
  install_dir: sage_install_dir / 'geometry',
)
subdir('polyhedron')
install_subdir(
  'riemannian_manifolds',
  install_dir: sage_install_dir / 'geometry',
)
subdir('triangulation')
