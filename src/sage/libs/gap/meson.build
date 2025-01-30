py.install_sources(
  '__init__.py',
  'all.py',
  'all_documented_functions.py',
  'assigned_names.py',
  'context_managers.py',
  'element.pxd',
  'gap_functions.py',
  'gap_globals.py',
  'gap_includes.pxd',
  'operations.py',
  'sage.gaprc',
  'saved_workspace.py',
  'test.py',
  'test_long.py',
  'util.pxd',
  subdir: 'sage/libs/gap',
)

# Ensure that the gaprc file is installed also in editable mode
fs.copyfile('sage.gaprc')

extension_data = {
  'element' : files('element.pyx'),
  'libgap' : files('libgap.pyx'),
  'util' : files('util.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gap, gmp]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/gap',
    install: true,
    include_directories: [inc_cpython, inc_rings],
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
      install_dir: sage_install_dir / 'libs/gap',
      install: true,
      copy: true,
    )
  endif
endforeach

