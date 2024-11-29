py.install_sources(
  'all__sagemath_objects.py',
  'ccobject.h',
  'cplusplus.pxd',
  'fast_callable.pxd',
  'fast_eval.pxd',
  'mod_int.h',
  'mod_int.pxd',
  'stdsage.pxd',
  subdir: 'sage/ext',
)

extension_data = {
  'fast_callable' : files('fast_callable.pyx'),
  'fast_eval' : files('fast_eval.pyx'),
  'memory' : files('memory.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/ext',
    install: true,
    include_directories: [inc_cpython, inc_ext],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

subdir('interpreters')
