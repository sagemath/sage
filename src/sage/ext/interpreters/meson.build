
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
    'wrapper_cdf.pyx',
    'wrapper_el.pyx',
    'wrapper_py.pyx',
    'wrapper_rdf.pyx',
    'wrapper_rr.pyx',
    '__init__.py',
  ],
  input: '../../../sage_setup/autogen/interpreters/internal/__init__.py',
  command: ['../../../sage_setup/autogen/interpreters/__main__.py', '@OUTDIR@'],
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
  'wrapper_cdf' : interpreters[8],
  'wrapper_el' : interpreters[9],
  'wrapper_py' : interpreters[10],
  'wrapper_rdf' : interpreters[11],
  'wrapper_rr' : interpreters[12],
}

# Manually create header files, which otherwise are not found
wrapper_el_header = custom_target(
  'wrapper_el.h',
  output: 'wrapper_el.h',
  input: [interpreters[9], interpreters],
  command: [
    cython.cmd_array(),
    '@INPUT0@',
    '-o',
    '@OUTPUT@',
    '-I',
    join_paths(meson.current_source_dir(), '../../../'),
  ],
)
wrapper_cc_header = custom_target(
  'wrapper_cc.h',
  output: 'wrapper_cc.h',
  input: [interpreters[7], interpreters],
  command: [
    cython.cmd_array(),
    '@INPUT0@',
    '-o',
    '@OUTPUT@',
    '-I',
    join_paths(meson.current_source_dir(), '../../../'),
  ],
)
wrapper_cdf_header = custom_target(
  'wrapper_cdf.h',
  output: 'wrapper_cdf.h',
  input: [interpreters[8], interpreters],
  command: [
    cython.cmd_array(),
    '@INPUT0@',
    '-o',
    '@OUTPUT@',
    '-I',
    join_paths(meson.current_source_dir(), '../../../'),
  ],
)
wrapper_rr_header = custom_target(
  'wrapper_rr.h',
  output: 'wrapper_rr.h',
  input: [interpreters[12], interpreters],
  command: [
    cython.cmd_array(),
    '@INPUT0@',
    '-o',
    '@OUTPUT@',
    '-I',
    join_paths(meson.current_source_dir(), '../../../'),
  ],
)

interpreters_dep = declare_dependency(
  sources: [
    wrapper_el_header,
    wrapper_cc_header,
    wrapper_cdf_header,
    wrapper_rr_header,
  ],
  include_directories: include_directories('.'),
)


foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/ext/interpreters',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_rings, inc_src],
    dependencies: [
      py_dep,
      cypari2,
      cysignals,
      gmp,
      gsl,
      mpc,
      mpfr,
      pari,
      interpreters_dep,
    ],
  )
endforeach
