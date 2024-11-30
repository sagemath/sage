brial = cc.find_library('brial', required: false, disabler: true)
# Cannot be found via pkg-config
brial_groebner = cc.find_library('brial_groebner', required: false, disabler: true)

py.install_sources(
  'PyPolyBoRi.py',
  '__init__.py',
  'blocks.py',
  'cnf.py',
  'easy_polynomials.py',
  'fglm.py',
  'frontend.py',
  'gbcore.py',
  'gbrefs.py',
  'heuristics.py',
  'interpolate.py',
  'interred.py',
  'll.py',
  'nf.py',
  'parallel.py',
  'pbori.pxd',
  'randompoly.py',
  'rank.py',
  'specialsets.py',
  'statistics.py',
  subdir: 'sage/rings/polynomial/pbori',
)

extension_data_cpp = {'pbori': files('pbori.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/polynomial/pbori',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_ext, inc_rings, inc_src],
    dependencies: [py_dep, brial, brial_groebner, cysignals, gmp, m4ri, png],
  )
endforeach

