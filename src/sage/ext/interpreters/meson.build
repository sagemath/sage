py.install_sources(
    'all.py',
    subdir: 'sage/ext/interpreters',
)

extension_data = {
    'wrapper_cc': files('wrapper_cc.pyx'),
    'wrapper_cdf': files('wrapper_cdf.pyx'),
    'wrapper_el': files('wrapper_el.pyx'),
    'wrapper_py': files('wrapper_py.pyx'),
    'wrapper_rdf': files('wrapper_rdf.pyx'),
    'wrapper_rr': files('wrapper_rr.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/ext/interpreters',
        install: true,
        include_directories: [inc_cpython, inc_ext, inc_interpreters, inc_rings, inc_src],
        dependencies: [py_dep, cypari2, cysignals, gmp, gsl, mpc, mpfr, pari, interpreters_dep],
    )
endforeach

