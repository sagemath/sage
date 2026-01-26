coxeter3 = cc.find_library(
  'coxeter3',
  required: get_option('coxeter3'),
  disabler: true,
)
py.install_sources(
  '__init__.py',
  'coxeter.pxd',
  'coxeter.pyx',
  'coxeter_group.py',
  'decl.pxd',
  subdir: 'sage/libs/coxeter3',
)

extension_data_cpp = {'coxeter': files('coxeter.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/coxeter3',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython],
    dependencies: [py_dep, cysignals, coxeter3],
  )
endforeach
