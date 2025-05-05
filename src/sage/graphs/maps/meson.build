py.install_sources(
  '__init__.py',
  'all.py',
  'banner.py',
  'custom_swap.py',
  'cycle_utils_provider.py',
  'cyclic_chained_list.py',
  'dynamic_planar_map_show.py',
  'example.py',
  'labelled_map.py',
  'map_decorator.py',
  'map_error.py',
  'map_generator.py',
  'map_permutation.py',
  'mutable_labelled_map.py',
  'mutable_topological_demi_edge.py',
  'permutation_utils_abstractor.py',
  'primitive_mutable_topological_demi_edge.py',
  'primitive_rotating_permutation.py',
  'primitive_rotating_permutation_utils_abstractor.py',
  'rooted_map.py',
  'rotating_permutation.py',
  'rotating_permutation_utils_abstractor.py',
  'splay_tree.py',
  'topological_demi_edge.py',
  subdir: 'sage/graphs/maps',
)

extension_data = {}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs/maps',
    install: true,
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

