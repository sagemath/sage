py.install_sources(
  'all.py',
  'automatic_semigroup.py',
  'free_abelian_monoid.py',
  'free_abelian_monoid_element.pxd',
  'free_monoid.py',
  'free_monoid_element.py',
  'hecke_monoid.py',
  'indexed_free_monoid.py',
  'monoid.py',
  'string_monoid.py',
  'string_monoid_element.py',
  'string_ops.py',
  'trace_monoid.py',
  subdir: 'sage/monoids',
)

extension_data = {
  'free_abelian_monoid_element' : files('free_abelian_monoid_element.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/monoids',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

