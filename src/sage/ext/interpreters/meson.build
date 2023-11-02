py.install_sources(
    'all.py',
    subdir: 'sage/ext/interpreters',
)

extension_data = {
    #'wrapper_cc': files('wrapper_cc.pyx'), # compilation error
    #'wrapper_cdf': files('wrapper_cdf.pyx'), # compilation error
    #'wrapper_el': files('wrapper_el.pyx'), # compilation error
    'wrapper_py': files('wrapper_py.pyx'),
    'wrapper_rdf': files('wrapper_rdf.pyx'),
    #'wrapper_rr': files('wrapper_rr.pyx'), # compilation error
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/ext/interpreters',
        install: true,
        include_directories: [inc_cpython, inc_ext, inc_interpreters, inc_rings, inc_src],
        dependencies: [py_dep, cypari2, cysignals, gmp, gsl, mpc, mpfr, pari],
    )
endforeach

