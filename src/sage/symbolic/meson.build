inc_ginac = include_directories('ginac')
inc_pynac = include_directories('.')

py.install_sources(
  'all.py',
  'assumptions.py',
  'benchmark.py',
  'callable.py',
  'complexity_measures.py',
  'constants.py',
  'expression.pxd',
  'expression_conversion_algebraic.py',
  'expression_conversion_sympy.py',
  'expression_conversions.py',
  'function.pxd',
  'function_factory.py',
  'maxima_wrapper.py',
  'operators.py',
  'pynac_wrap.h',
  'random_tests.py',
  'relation.py',
  'ring.pxd',
  'subring.py',
  'symbols.py',
  'symengine.py',
  'tests.py',
  'units.py',
  subdir: 'sage/symbolic',
)

extension_data = {
  'function' : files('function.pyx'),
  'ring' : files('ring.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/symbolic',
    install: true,
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_ginac,
      inc_gsl,
      inc_pynac,
      inc_rings,
    ],
    dependencies: [
      py_dep,
      cysignals,
      gmp,
    ],
  )
endforeach

extension_data_cpp = {
  'expression': files(
    'expression.pyx',
    'ginac/add.cpp',
    'ginac/archive.cpp',
    'ginac/assume.cpp',
    'ginac/basic.cpp',
    'ginac/cmatcher.cpp',
    'ginac/constant.cpp',
    'ginac/context.cpp',
    'ginac/ex.cpp',
    'ginac/expair.cpp',
    'ginac/expairseq.cpp',
    'ginac/exprseq.cpp',
    'ginac/fderivative.cpp',
    'ginac/function.cpp',
    'ginac/function_info.cpp',
    'ginac/infinity.cpp',
    'ginac/infoflagbase.cpp',
    'ginac/inifcns.cpp',
    'ginac/inifcns_comb.cpp',
    'ginac/inifcns_gamma.cpp',
    'ginac/inifcns_hyperb.cpp',
    'ginac/inifcns_hyperg.cpp',
    'ginac/inifcns_nstdsums.cpp',
    'ginac/inifcns_orthopoly.cpp',
    'ginac/inifcns_trans.cpp',
    'ginac/inifcns_trig.cpp',
    'ginac/inifcns_zeta.cpp',
    'ginac/lst.cpp',
    'ginac/matrix.cpp',
    'ginac/mpoly-ginac.cpp',
    'ginac/mpoly-singular.cpp',
    'ginac/mpoly.cpp',
    'ginac/mul.cpp',
    'ginac/normal.cpp',
    'ginac/numeric.cpp',
    'ginac/operators.cpp',
    'ginac/order.cpp',
    'ginac/power.cpp',
    'ginac/print.cpp',
    'ginac/pseries.cpp',
    'ginac/py_funcs.cpp',
    'ginac/registrar.cpp',
    'ginac/relational.cpp',
    'ginac/remember.cpp',
    'ginac/sum.cpp',
    'ginac/symbol.cpp',
    'ginac/templates.cpp',
    'ginac/upoly-ginac.cpp',
    'ginac/useries.cpp',
    'ginac/utils.cpp',
    'ginac/wildcard.cpp',
  ),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/symbolic',
    install: true,
    override_options: ['cython_language=cpp'],
    cpp_args: '-std=c++11',
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_ginac,
      inc_gsl,
      inc_pynac,
      inc_rings,
      include_directories('../libs/gmp'),
    ],
    dependencies: [
      py_dep,
      cysignals,
      gmp,
      gsl,
      singular_factory,
    ],
  )
endforeach

install_subdir('integration', install_dir: sage_install_dir / 'symbolic')
