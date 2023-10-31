py.install_sources(
    'all.py',
    subdir: 'sage/ext/interpreters',
)

extension_data = {
    #'wrapper_cc': files('wrapper_cc.pyx'), - compiliation error
    #'wrapper_cdf': files('wrapper_cdf.pyx'), - compiliation error
    #'wrapper_el': files('wrapper_el.pyx'), - compiliation error
    'wrapper_py': files('wrapper_py.pyx'),
    'wrapper_rdf': files('wrapper_rdf.pyx'),
    #'wrapper_rr': files('wrapper_rr.pyx'), - compiliation error
}
extension_data_cpp = {
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/ext/interpreters',
        install: true,
        include_directories: [inc_cpython, inc_src],
        dependencies: [py_dep, cypari2, cysignals],
    )
endforeach
foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/ext/interpreters',
        install: true,
        override_options: ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_src],
        dependencies: [py_dep, cypari2, cysignals],
    )
endforeach
