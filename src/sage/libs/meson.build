sirocco = cc.find_library('sirocco', required: false, disabler: true)
# cannot be found via pkg-config
ecl = cc.find_library('ecl')
braiding = cc.find_library('braiding')
gc = cc.find_library('gc')
homfly = cc.find_library('homfly', has_headers: ['homfly.h'])

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

dependencies = [py_dep, braiding, cysignals, ecl, ecm, gc, gmp, homfly]

foreach name, pyx : extension_data
  deps = dependencies
  if name == 'sirocco'
    deps += [sirocco]
  elif name == 'meataxe'
    deps += [mtx]
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
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/libs',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
    dependencies: dependencies,
  )
endforeach

subdir('arb')
subdir('coxeter3')
install_subdir('cremona', install_dir: sage_install_dir / 'libs')
subdir('eclib')
subdir('flint')
subdir('gap')
subdir('giac')
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
