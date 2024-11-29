py.install_sources(
  '__init__.py',
  'acb.pxd',
  'acb_calc.pxd',
  'acb_elliptic.pxd',
  'acb_hypgeom.pxd',
  'acb_mat.pxd',
  'acb_modular.pxd',
  'acb_poly.pxd',
  'arb.pxd',
  'arb_fmpz_poly.pxd',
  'arb_hypgeom.pxd',
  'arb_wrap.h',
  'arf.pxd',
  'bernoulli.pxd',
  'mag.pxd',
  'types.pxd',
  subdir: 'sage/libs/arb',
)

extension_data = {'arith' : files('arith.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/arb',
    install: true,
    include_directories: [inc_cpython, inc_flint, inc_rings],
    dependencies: [py_dep, flint, gmp, mpfr],
  )
endforeach

