py.install_sources('all.py', 'random_variable.py', subdir: 'sage/probability')

extension_data = {
  'probability_distribution' : files('probability_distribution.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/probability',
    install: true,
    include_directories: [],
    dependencies: [py_dep, cysignals, gmp, gsl],
  )
endforeach

