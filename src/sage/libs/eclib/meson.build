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
  deps = [py_dep, cysignals, ec, flint, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/eclib',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ext, inc_flint, inc_ntl, inc_rings],
    dependencies: deps,
  )

  all_deps_found = true
  foreach dep : deps
    if not dep.found()
      all_deps_found = false
      break
    endif
  endforeach
  if not all_deps_found
    configure_file(
      input: not_available,
      output: name + '.py',
      install_tag: 'python-runtime',
      install_dir: sage_install_dir / 'libs/eclib',
      install: true,
      copy: true,
    )
  endif
endforeach

