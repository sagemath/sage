py.install_sources(
    'all.py',
    'all__sagemath_coxeter3.py',
    'all__sagemath_meataxe.py',
    'all__sagemath_objects.py',
    'all__sagemath_sirocco.py',
    subdir: 'sage/libs',
)

extension_data_cpp = {
    'braiding': files('braiding.pyx'),
}
extension_data = {
    'ecl': files('ecl.pyx'),
    'homfly': files('homfly.pyx'),
    'libecm': files('libecm.pyx'),
    #'meataxe': files('meataxe.pyx'),# - need meataxe
}

dependencies = [py_dep, cysignals, cypari2]

sirocco = dependency('sirocco', required: false)
if sirocco.found()
    #extension_data['sirocco'] += files('sirocco.pxd')
    dependencies += [sirocco]
endif

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs',
        install: true,
        include_directories: [inc_cpython],
        dependencies: dependencies,
    )
endforeach
foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython],
        dependencies: dependencies,
    )
endforeach

install_subdir('cremona', install_dir: sage_install_dir / 'libs')
install_subdir('linkages', install_dir: sage_install_dir / 'libs')
install_subdir('lrcalc', install_dir: sage_install_dir / 'libs')
install_subdir('mwrank', install_dir: sage_install_dir / 'libs')
install_subdir('pynac', install_dir: sage_install_dir / 'libs')

subdir('arb')
subdir('coxeter3')
subdir('eclib')
subdir('flint')
subdir('gap')
subdir('giac')
subdir('glpk')
subdir('gmp')
subdir('gsl')
subdir('lcalc')
subdir('linbox')
subdir('mpmath')
subdir('ntl')
subdir('pari')
subdir('singular')
subdir('symmetrica')
