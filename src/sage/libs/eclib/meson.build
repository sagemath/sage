py.install_sources(
  '__init__.pxd',
  '__init__.py',
  'all.py',
  'constructor.py',
  'homspace.pxd',
  'interface.py',
  'mat.pxd',
  'newforms.pxd',
  subdir: 'sage/libs/eclib',
)

extension_data_cpp = {
  'homspace': files('homspace.pyx'),
  'mat': files('mat.pyx'),
  'mwrank': files('mwrank.pyx'),
  'newforms': files('newforms.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/eclib',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_ntl, inc_rings],
    dependencies: [py_dep, cysignals, ec, ecl, flint, gmp],
  )
endforeach

