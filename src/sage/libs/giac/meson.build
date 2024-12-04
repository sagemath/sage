giac = cc.find_library('giac', required: false, disabler: true)

py.install_sources('__init__.py', 'giac.pxd', 'misc.h', subdir: 'sage/libs/giac')

extension_data_cpp = {'giac': files('giac.pyx')}

cpp_args = []
if is_windows
  # Giac doesn't compile with c++17 on Windows: https://github.com/sagemath/sage/issues/38985
  cpp_args += '-std=c++11'
endif

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs/giac',
    install: true,
    override_options: ['cython_language=cpp'],
    cpp_args: cpp_args,
    include_directories: [inc_cpython, inc_ext, inc_rings],
    dependencies: [py_dep, cysignals, giac, gmp],
  )
endforeach
