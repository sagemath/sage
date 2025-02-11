py.install_sources(
  '__init__.py',
  'all.py',
  'base.pxd',
  'implicit_plot3d.py',
  'index_face_set.pxd',
  'introduction.py',
  'list_plot3d.py',
  'parametric_plot3d.py',
  'parametric_surface.pxd',
  'platonic.py',
  'plot3d.py',
  'plot_field3d.py',
  'revolution_plot3d.py',
  'shapes.pxd',
  'shapes2.py',
  'tachyon.py',
  'texture.py',
  'transform.pxd',
  'tri_plot.py',
  subdir: 'sage/plot/plot3d',
)

extension_data = {
  'base' : files('base.pyx'),
  'implicit_surface' : files('implicit_surface.pyx'),
  'index_face_set' : files('index_face_set.pyx'),
  'shapes' : files('shapes.pyx'),
  'transform' : files('transform.pyx'),
  'parametric_surface' : files('parametric_surface.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp]
  if name == 'parametric_surface'
    deps += [interpreters_dep]
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/plot/plot3d',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_numpy],
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
      install_dir: sage_install_dir / 'plot/plot3d',
      install: true,
      copy: true,
    )
  endif
endforeach
