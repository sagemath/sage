py.install_sources(
    'all.py',
    'f_matrix.py',
    'fusion_double.py',
    'fusion_ring.py',
    subdir: 'sage/algebras/fusion_rings',
)

extension_data = {
    'fast_parallel_fusion_ring_braid_repn': files('fast_parallel_fusion_ring_braid_repn.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/algebras/fusion_rings',
        install: true,
        include_directories: [inc_numpy, inc_ntl],
        dependencies: [py_dep, cysignals],
    )
endforeach

extension_data_cpp = {
    'fast_parallel_fmats_methods': files('fast_parallel_fmats_methods.pyx'),
    'shm_managers': files('shm_managers.pyx'),
    'poly_tup_engine': files('poly_tup_engine.pyx'),
}
foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/algebras/fusion_rings',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_numpy, inc_ntl],
        dependencies: [py_dep, cysignals],
    )
endforeach
