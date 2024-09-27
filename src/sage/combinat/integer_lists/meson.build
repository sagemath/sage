py.install_sources(
  '__init__.py',
  'lists.py',
  'nn.py',
  subdir: 'sage/combinat/integer_lists',
)

extension_data = {'base' : files('base.pyx'), 'invlex' : files('invlex.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/integer_lists',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, gmp],
  )
endforeach

