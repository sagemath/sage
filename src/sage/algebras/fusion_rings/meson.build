py.install_sources(
  '__init__.py',
  'all.py',
  'f_matrix.py',
  'fast_parallel_fmats_methods.pxd',
  'fast_parallel_fusion_ring_braid_repn.pxd',
  'fusion_double.py',
  'fusion_ring.py',
  'poly_tup_engine.pxd',
  'shm_managers.pxd',
  subdir: 'sage/algebras/fusion_rings',
)

extension_data = {
  'fast_parallel_fusion_ring_braid_repn' : files(
    'fast_parallel_fusion_ring_braid_repn.pyx',
  ),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/algebras/fusion_rings',
    install: true,
    include_directories: [inc_cpython, inc_ntl, inc_rings],
    dependencies: [py_dep, cysignals, gmp, numpy],
  )
endforeach

extension_data_cpp = {
  'fast_parallel_fmats_methods': files('fast_parallel_fmats_methods.pyx'),
  'poly_tup_engine': files('poly_tup_engine.pyx'),
  'shm_managers': files('shm_managers.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/algebras/fusion_rings',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ntl, inc_rings],
    dependencies: [py_dep, cysignals, gmp, numpy, singular],
  )
endforeach

