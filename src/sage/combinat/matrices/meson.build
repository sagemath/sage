py.install_sources(
  '__init__.py',
  'all.py',
  'dlxcpp.py',
  'hadamard_matrix.py',
  'latin.py',
  subdir: 'sage/combinat/matrices',
)

extension_data_cpp = {'dancing_links': files('dancing_links.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/matrices',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

