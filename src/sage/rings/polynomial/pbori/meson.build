brial = dependency('brial', required: false, disabler: true)
if not brial.found()
  # fallback for older versions without pkg-config
  brial = cc.find_library('brial', required: get_option('brial'), disabler: true)
endif
# Cannot be found via pkg-config
brial_groebner = cc.find_library(
  'brial_groebner',
  required: get_option('brial'),
  disabler: true,
)

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
  'pbori.pyx',
  'randompoly.py',
  'rank.py',
  'specialsets.py',
  'statistics.py',
  subdir: 'sage/rings/polynomial/pbori',
)

extension_data_cpp = {'pbori': files('pbori.pyx')}

foreach name, pyx : extension_data_cpp
  override_options = ['cython_language=cpp']
  if not boost_cpp17_compatible
    override_options += ['cpp_std=c++11']
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/polynomial/pbori',
    install: true,
    override_options: override_options,
    include_directories: [inc_cpython, inc_ext, inc_rings, inc_src],
    dependencies: [
      py_dep,
      boost,
      brial,
      brial_groebner,
      cysignals,
      gmp,
      m4ri,
      png,
    ],
  )
endforeach

