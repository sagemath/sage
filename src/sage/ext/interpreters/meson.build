
interpreters = custom_target(
  'sage.ext.interpreters',
  output: [
    'all.py',
    'wrapper_cc.pxd',
    'wrapper_cdf.pxd',
    'wrapper_el.pxd',
    'wrapper_py.pxd',
    'wrapper_rdf.pxd',
    'wrapper_rr.pxd',
    'wrapper_cc.pyx',
    'interp_cc.c',
    'wrapper_cdf.pyx',
    'interp_cdf.c',
    'wrapper_el.pyx',
    'interp_el.c',
    'wrapper_py.pyx',
    'interp_py.c',
    'wrapper_rdf.pyx',
    'interp_rdf.c',
    'wrapper_rr.pyx',
    'interp_rr.c',
    '__init__.py',
  ],
  input: '../../../autogen/interpreters/internal/__init__.py',
  command: ['../../../autogen/interpreters/__main__.py', '@OUTDIR@'],
  # Manually install the generated files instead of using install_sources
  # this is a workaround for https://github.com/mesonbuild/meson/issues/7372
  install: true,
  install_dir: py.get_install_dir() / 'sage/ext/interpreters',
  install_tag: 'python-runtime',
)

# Use this once https://github.com/mesonbuild/meson/issues/7372 is fixed
#foreach file : interpreters.to_list()
#  py.install_sources(
#    file.full_path(),
#    subdir : 'sage/ext/interpreters'
#  )
#endforeach

extension_data = {
  'wrapper_cc' : interpreters[7],
  'wrapper_cdf' : interpreters[9],
  'wrapper_el' : interpreters[11],
  'wrapper_py' : interpreters[13],
  'wrapper_rdf' : interpreters[15],
  'wrapper_rr' : interpreters[17],
}

interpreters_dep = declare_dependency(
  include_directories: include_directories('.'),
  sources: [
    interpreters[1],
    interpreters[2],
    interpreters[3],
    interpreters[4],
    interpreters[5],
  ],
)

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/ext/interpreters',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_rings, inc_src],
    dependencies: [py_dep, cypari2, cysignals, gmp, gsl, mpc, mpfr, pari],
  )
endforeach
