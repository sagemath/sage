py.install_sources(
  '__init__.py',
  'all.py',
  'arithgroup_generic.py',
  'arithgroup_perm.py',
  'congroup_gamma.py',
  'congroup_gamma0.py',
  'congroup_gamma1.py',
  'congroup_gammaH.py',
  'congroup_generic.py',
  'congroup_sl2z.py',
  'tests.py',
  subdir: 'sage/modular/arithgroup',
)

extension_data = {
  'arithgroup_element' : files('arithgroup_element.pyx'),
  'congroup' : files('congroup.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modular/arithgroup',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_rings],
    dependencies: [py_dep, cysignals, flint, gmp],
  )
endforeach

# Manually create header file, which otherwise is not found
farey_symbol_header = custom_target(
  'farey_symbol.h',
  output: 'farey_symbol.h',
  input: 'farey_symbol.pyx',
  command: [
    cython.cmd_array(),
    '--cplus',
    '@INPUT@',
    '-o',
    '@OUTPUT@',
    '-I',
    join_paths(meson.current_source_dir(), '../../../'),
  ],
)

extension_data_cpp = {
  'farey_symbol': [
    files('farey.cpp', 'farey_symbol.pyx', 'sl2z.cpp'),
    farey_symbol_header[0],
  ],
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modular/arithgroup',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_rings, inc_src],
    dependencies: [py_dep, cysignals, flint, gmp, gmpxx],
  )
endforeach

