py.install_sources(
  '__init__.py',
  'all.pxd',
  'binop.pxd',
  'misc.pxd',
  'mpf.pxd',
  'mpn.pxd',
  'mpq.pxd',
  'mpz.pxd',
  'pylong.pxd',
  'random.pxd',
  'randomize.pxd',
  'types.pxd',
  subdir: 'sage/libs/gmp',
)

extension_data = {'pylong' : files('pylong.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/gmp',
    install: true,
    include_directories: [inc_cpython],
    dependencies: [py_dep, gmp],
  )
endforeach

