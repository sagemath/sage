py.install_sources(
  '__init__.py',
  'all.py',
  'animate.py',
  'arc.py',
  'arrow.py',
  'bar_chart.py',
  'bezier_path.py',
  'circle.py',
  'colors.py',
  'contour_plot.py',
  'density_plot.py',
  'disk.py',
  'ellipse.py',
  'graphics.py',
  'histogram.py',
  'hyperbolic_arc.py',
  'hyperbolic_polygon.py',
  'hyperbolic_regular_polygon.py',
  'line.py',
  'matrix_plot.py',
  'misc.py',
  'multigraphics.py',
  'plot.py',
  'plot_field.py',
  'point.py',
  'polygon.py',
  'primitive.py',
  'scatter_plot.py',
  'step.py',
  'streamline_plot.py',
  'text.py',
  subdir: 'sage/plot',
)

extension_data = {'complex_plot' : files('complex_plot.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/plot',
    install: true,
    include_directories: [inc_cpython, inc_gsl, inc_rings],
    dependencies: [py_dep, cysignals, gmp, gsl, numpy],
  )
endforeach

subdir('plot3d')
