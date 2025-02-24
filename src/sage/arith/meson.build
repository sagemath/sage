py.install_sources(
  'all.py',
  'all__sagemath_objects.py',
  'constants.pxd',
  'functions.pxd',
  'long.pxd',
  'misc.py',
  'multi_modular.pxd',
  'numerical_approx.pxd',
  'power.pxd',
  'rational_reconstruction.pxd',
  subdir: 'sage/arith',
)

extension_data = {
  'functions' : files('functions.pyx'),
  'multi_modular' : files('multi_modular.pyx'),
  'numerical_approx' : files('numerical_approx.pyx'),
  'power' : files('power.pyx'),
  'rational_reconstruction' : files('rational_reconstruction.pyx'),
  'srange' : files('srange.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/arith',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

