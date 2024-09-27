py.install_sources(
  'all.py',
  'ambient.py',
  'boundary.py',
  'element.py',
  'g1list.py',
  'ghlist.py',
  'hecke_operator.py',
  'manin_symbol_list.py',
  'modsym.py',
  'modular_symbols.py',
  'p1list_nf.py',
  'relation_matrix.py',
  'space.py',
  'subspace.py',
  'tests.py',
  subdir: 'sage/modular/modsym',
)

extension_data = {
  'apply' : files('apply.pyx'),
  'heilbronn' : files('heilbronn.pyx'),
  'manin_symbol' : files('manin_symbol.pyx'),
  'p1list' : files('p1list.pyx'),
  'relation_matrix_pyx' : files('relation_matrix_pyx.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modular/modsym',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_rings],
    dependencies: [py_dep, cysignals, flint, gmp],
  )
endforeach

