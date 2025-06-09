sirocco = cc.find_library('sirocco', required: false, disabler: true)
# cannot be found via pkg-config
ecl = cc.find_library('ecl', required: false, disabler: true)
if not ecl.found() and not is_windows
  ecl_proj = subproject('ecl')
  ecl = ecl_proj.get_variable('ecl_dep')
endif
braiding = dependency('libbraiding', required: false)
if not braiding.found()
  # Fallback since pkg-config support was only added in v1.3.1
  braiding = cc.find_library(
    'braiding',
    required: not is_windows,
    disabler: true,
  )
endif
gc = dependency(
  ['bdw-gc-threaded', 'bdw-gc'],
  version: '>=7.6.4',
  required: not is_windows,
  disabler: true,
)
homfly = cc.find_library(
  'homfly',
  has_headers: ['homfly.h'],
  required: false,
  disabler: true,
)

py.install_sources(
  '__init__.py',
  'all.py',
  'all__sagemath_coxeter3.py',
  'all__sagemath_meataxe.py',
  'all__sagemath_objects.py',
  'all__sagemath_sirocco.py',
  'ecl.pxd',
  'eclsig.h',
  'gmpxx.pxd',
  'iml.pxd',
  'm4ri.pxd',
  'm4rie.pxd',
  'meataxe.pxd',
  subdir: 'sage/libs',
)

extension_data = {
  'ecl' : files('ecl.pyx'),
  'homfly' : files('homfly.pyx'),
  'libecm' : files('libecm.pyx'),
  'sirocco': files('sirocco.pyx'),
  'meataxe': files('meataxe.pyx'),
}

dependencies = [py_dep, cysignals, ecl, gc, gmp]

foreach name, pyx : extension_data
  deps = dependencies
  if name == 'sirocco'
    deps += [sirocco]
  elif name == 'meataxe'
    deps += [mtx]
  elif name == 'homfly'
    deps += [homfly]
  elif name == 'libecm'
    deps += [ecm]
  endif

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: deps,
  )
endforeach

extension_data_cpp = {'braiding': files('braiding.pyx')}

foreach name, pyx : extension_data_cpp
  deps = dependencies
  if name == 'braiding'
    deps += [braiding]
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
    dependencies: deps,
  )
endforeach

subdir('arb')
subdir('coxeter3')
install_subdir('cremona', install_dir: sage_install_dir / 'libs')
subdir('eclib')
subdir('flint')
subdir('gap')
install_subdir('giac', install_dir: sage_install_dir / 'libs')
subdir('glpk')
subdir('gmp')
subdir('gsl')
subdir('lcalc')
subdir('linbox')
install_subdir('linkages', install_dir: sage_install_dir / 'libs')
install_subdir('lrcalc', install_dir: sage_install_dir / 'libs')
install_subdir('mpfr', install_dir: sage_install_dir / 'libs')
subdir('mpmath')
install_subdir('mwrank', install_dir: sage_install_dir / 'libs')
subdir('ntl')
subdir('pari')
subdir('singular')
subdir('symmetrica')
