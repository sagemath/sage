tdlib = dependency('tdlib', required: false)

py.install_sources(
    'all.py',
    'all__sagemath_tdlib.py',
    'modular_decomposition.py',
    subdir: 'sage/graphs/graph_decompositions',
)

extension_data = {
    'bandwidth': files('bandwidth.pyx'),
    'cutwidth': files('cutwidth.pyx'),
    'fast_digraph': files('fast_digraph.pyx'),
    'graph_products': files('graph_products.pyx'),
    'rankwidth': files('rankwidth.pyx'),
    'tree_decomposition': files('tree_decomposition.pyx'),
    'vertex_separation': files('vertex_separation.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/graphs/graph_decompositions',
        install: true,
        include_directories: [inc_cpython, inc_data_structures],
        dependencies: [py_dep, cysignals],
    )
endforeach

extension_data_cpp = {
    'clique_separators': files('clique_separators.pyx'),
}

foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/graphs/graph_decompositions',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_data_structures],
        dependencies: [py_dep, cysignals],
    )
endforeach

if tdlib.found()
    py.extension_module('tdlib',
        sources: 'tdlib.pyx',
        subdir: 'sage/graphs/graph_decompositions',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_data_structures],
        dependencies: [py_dep, cysignals, tdlib],
    )
endif

