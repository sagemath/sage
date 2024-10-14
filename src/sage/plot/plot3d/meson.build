py.install_sources(
  'all.py',
  'implicit_plot3d.py',
  'introduction.py',
  'list_plot3d.py',
  'parametric_plot3d.py',
  'platonic.py',
  'plot3d.py',
  'plot_field3d.py',
  'revolution_plot3d.py',
  'shapes2.py',
  'tachyon.py',
  'texture.py',
  'tri_plot.py',
  subdir: 'sage/plot/plot3d',
)

extension_data = {
  'base' : files('base.pyx'),
  'implicit_surface' : files('implicit_surface.pyx'),
  'index_face_set' : files('index_face_set.pyx'),
  'shapes' : files('shapes.pyx'),
  'transform' : files('transform.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/plot/plot3d',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_numpy],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

# This extension needs the interpreters, but the others don't.
py.extension_module(
  'parametric_surface',
  sources: files('parametric_surface.pyx'),
  subdir: 'sage/plot/plot3d',
  install: true,
  include_directories: [inc_cpython, inc_ext, inc_numpy],
  dependencies: [py_dep, cysignals, gmp, interpreters_dep],
)
