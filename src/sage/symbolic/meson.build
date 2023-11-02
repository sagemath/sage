# cannot be found via pkg-config
singular = declare_dependency(
  dependencies: [
    cc.find_library('singular', required: false),
  ],
)
ginac = declare_dependency(
  dependencies: [
    cc.find_library('ginac', required: false),
  ],
)
inc_ginac = include_directories('ginac')
inc_pynac = include_directories('.')

py.install_sources(
    '__init__.py',
    'all.py',
    'assumptions.py',
    'benchmark.py',
    'callable.py',
    'complexity_measures.py',
    'constants.py',
    'expression_conversion_algebraic.py',
    'expression_conversion_sympy.py',
    'expression_conversions.py',
    'function_factory.py',
    'maxima_wrapper.py',
    'operators.py',
    'random_tests.py',
    'relation.py',
    'subring.py',
    'symbols.py',
    'symengine.py',
    'tests.py',
    'units.py',
    subdir: 'sage/symbolic',
)

extension_data = {
    'function': files('function.pyx'),
    'ring': files('ring.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/symbolic',
        install: true,
        include_directories: [inc_cpython, inc_ext, inc_ginac, inc_gsl, inc_pynac, inc_rings],
        dependencies: [py_dep, cypari2, cysignals, flint, ginac, gmp, gsl, mpfr, pari, singular],
    )
endforeach

extension_data_cpp = {
    'expression': files('expression.pyx'),
}

foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/symbolic',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_ext, inc_ginac, inc_gsl, inc_pynac, inc_rings],
        dependencies: [py_dep, cypari2, cysignals, flint, ginac, gmp, gsl, mpfr, pari, singular],
    )
endforeach

install_subdir('integration', install_dir: sage_install_dir / 'symbolic')
