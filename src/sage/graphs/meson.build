bliss = cc.find_library('bliss', required: false, disabler: true)
mcqd = cc.find_library('mcqd', required: false, disabler: true)
cliquer = cc.find_library('cliquer')

# Cannot be found via pkg-config
planarity = cc.find_library('planarity')

py.install_sources(
  'all.py',
  'all__sagemath_bliss.py',
  'all__sagemath_mcqd.py',
  'all__sagemath_tdlib.py',
  'bipartite_graph.py',
  'cliquer.pxd',
  'cographs.py',
  'connectivity.pxd',
  'convexity_properties.pxd',
  'digraph.py',
  'digraph_generators.py',
  'distances_all_pairs.pxd',
  'domination.py',
  'dot2tex_utils.py',
  'generic_graph.py',
  'generic_graph_pyx.pxd',
  'graph.py',
  'graph_database.py',
  'graph_editor.py',
  'graph_generators.py',
  'graph_input.py',
  'graph_latex.py',
  'graph_list.py',
  'graph_plot.py',
  'graph_plot_js.py',
  'hypergraph_generators.py',
  'independent_sets.pxd',
  'isgci.py',
  'lovasz_theta.py',
  'matching.py',
  'matching_covered_graph.py',
  'mcqd.pxd',
  'orientations.py',
  'partial_cube.py',
  'pq_trees.py',
  'print_graphs.py',
  'schnyder.py',
  'traversals.pxd',
  'trees.pxd',
  'tutte_polynomial.py',
  subdir: 'sage/graphs',
)

extension_data = {
  'asteroidal_triples' : files('asteroidal_triples.pyx'),
  'centrality' : files('centrality.pyx'),
  'chrompoly' : files('chrompoly.pyx'),
  'cliquer' : files('cliquer.pyx'),
  'comparability' : files('comparability.pyx'),
  'connectivity' : files('connectivity.pyx'),
  'convexity_properties' : files('convexity_properties.pyx'),
  'distances_all_pairs' : files('distances_all_pairs.pyx'),
  'generic_graph_pyx' : files('generic_graph_pyx.pyx'),
  'genus' : files('genus.pyx'),
  'graph_generators_pyx' : files('graph_generators_pyx.pyx'),
  'hyperbolicity' : files('hyperbolicity.pyx'),
  'independent_sets' : files('independent_sets.pyx'),
  'isoperimetric_inequalities' : files('isoperimetric_inequalities.pyx'),
  'line_graph' : files('line_graph.pyx'),
  'matchpoly' : files('matchpoly.pyx'),
  'planarity' : files('planarity.pyx'),
  'spanning_tree' : files('spanning_tree.pyx'),
  'strongly_regular_db' : files('strongly_regular_db.pyx'),
  'trees' : files('trees.pyx'),
  'views' : files('views.pyx'),
  'weakly_chordal' : files('weakly_chordal.pyx'),
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
