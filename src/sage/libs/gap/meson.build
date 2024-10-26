py.install_sources(
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
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/gap',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cysignals, gap, gmp],
  )
endforeach

