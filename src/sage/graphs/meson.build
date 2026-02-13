bliss = cc.find_library('bliss', required: get_option('bliss'), disabler: true)

mcqd = disabler()
if cc.has_header('mcqd.h', required: get_option('mcqd'))
  # mcqd is a header-only library
  mcqd = declare_dependency()
endif

cliquer = dependency('libcliquer', required: false, disabler: true)
if not cliquer.found()
  # fallback for older versions without pkg-config
  cliquer = cc.find_library('cliquer', required: not is_windows, disabler: true)
endif

planarity = dependency('libplanarity', required: false, disabler: true)
if not planarity.found()
  # fallback for older versions without pkg-config
  planarity = cc.find_library(
    'planarity',
    required: not is_windows,
    disabler: true,
  )
endif

py.install_sources(
  '__init__.py',
  'all.py',
  'asteroidal_triples.pyx',
  'bipartite_graph.py',
  'bliss.pyx',
  'centrality.pyx',
  'chrompoly.pyx',
  'cliquer.pxd',
  'cliquer.pyx',
  'cographs.py',
  'comparability.pyx',
  'connectivity.pxd',
  'connectivity.pyx',
  'convexity_properties.pxd',
  'convexity_properties.pyx',
  'cycle_enumeration.py',
  'digraph.py',
  'digraph_generators.py',
  'distances_all_pairs.pxd',
  'distances_all_pairs.pyx',
  'domination.py',
  'dot2tex_utils.py',
  'edge_connectivity.pyx',
  'generic_graph.py',
  'generic_graph_pyx.pxd',
  'generic_graph_pyx.pyx',
  'genus.pyx',
  'graph.py',
  'graph_coloring.pyx',
  'graph_database.py',
  'graph_editor.py',
  'graph_generators.py',
  'graph_generators_pyx.pyx',
  'graph_input.py',
  'graph_latex.py',
  'graph_list.py',
  'graph_plot.py',
  'graph_plot_js.py',
  'hyperbolicity.pyx',
  'hypergraph_generators.py',
  'independent_sets.pxd',
  'independent_sets.pyx',
  'isgci.py',
  'isoperimetric_inequalities.pyx',
  'line_graph.pyx',
  'lovasz_theta.py',
  'matching.py',
  'matching_covered_graph.py',
  'matchpoly.pyx',
  'mcqd.pxd',
  'mcqd.pyx',
  'morphisms.py',
  'orientations.py',
  'partial_cube.py',
  'path_enumeration.pyx',
  'planarity.pyx',
  'pq_trees.py',
  'print_graphs.py',
  'schnyder.py',
  'spanning_tree.pyx',
  'strongly_regular_db.pyx',
  'traversals.pxd',
  'traversals.pyx',
  'tutte_polynomial.py',
  'views.pyx',
  'weakly_chordal.pyx',
  subdir: 'sage/graphs',
)


extension_data = {
  'asteroidal_triples': files('asteroidal_triples.pyx'),
  'centrality': files('centrality.pyx'),
  'chrompoly': files('chrompoly.pyx'),
  'cliquer': files('cliquer.pyx'),
  'comparability': files('comparability.pyx'),
  'connectivity': files('connectivity.pyx'),
  'convexity_properties': files('convexity_properties.pyx'),
  'distances_all_pairs': files('distances_all_pairs.pyx'),
  'generic_graph_pyx': files('generic_graph_pyx.pyx'),
  'genus': files('genus.pyx'),
  'graph_generators_pyx': files('graph_generators_pyx.pyx'),
  'hyperbolicity': files('hyperbolicity.pyx'),
  'independent_sets': files('independent_sets.pyx'),
  'isoperimetric_inequalities': files('isoperimetric_inequalities.pyx'),
  'line_graph': files('line_graph.pyx'),
  'matchpoly': files('matchpoly.pyx'),
  'planarity': files('planarity.pyx'),
  'spanning_tree': files('spanning_tree.pyx'),
  'strongly_regular_db': files('strongly_regular_db.pyx'),
  'views': files('views.pyx'),
  'weakly_chordal': files('weakly_chordal.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs',
    install: true,
    include_directories: [
      inc_cpython,
      inc_data_structures,
      inc_flint,
      inc_rings,
    ],
    dependencies: [py_dep, cysignals, cliquer, flint, gmp, planarity],
  )
endforeach

extension_data_cpp = {
  'edge_connectivity': files('edge_connectivity.pyx'),
  'graph_coloring': files('graph_coloring.pyx'),
  'path_enumeration': files('path_enumeration.pyx'),
  'traversals': files('traversals.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/graphs',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [
      inc_cpython,
      inc_data_structures,
      inc_flint,
      inc_rings,
    ],
    dependencies: [py_dep, cysignals, flint, gmp, planarity],
  )
endforeach

py.extension_module(
  'bliss',
  sources: files('bliss.pyx'),
  subdir: 'sage/graphs',
  install: true,
  override_options: ['cython_language=cpp'],
  include_directories: [inc_cpython, inc_data_structures, inc_flint, inc_rings],
  dependencies: [py_dep, cysignals, bliss],
)

py.extension_module(
  'mcqd',
  sources: files('mcqd.pyx'),
  subdir: 'sage/graphs',
  install: true,
  override_options: ['cython_language=cpp'],
  include_directories: [inc_cpython, inc_data_structures, inc_flint, inc_rings],
  dependencies: [py_dep, cysignals, mcqd],
)

subdir('base')
subdir('generators')
subdir('graph_decompositions')
subdir('tests')
