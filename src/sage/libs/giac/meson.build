giac = cc.find_library('giac', required: false, disabler: true)

py.install_sources('__init__.py', 'giac.pxd', 'misc.h', subdir: 'sage/libs/giac')

extension_data_cpp = {'giac': files('giac.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/giac',
    install: true,
    override_options: [
      'cython_language=cpp',
      # Temporary workaround for https://github.com/sagemath/sage/issues/38985
      'cpp_std=c++11',
    ],
    include_directories: [inc_cpython, inc_ext, inc_rings],
    dependencies: [py_dep, cysignals, giac, gmp],
  )
endforeach
