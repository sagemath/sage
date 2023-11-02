py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/lcalc',
)


extension_data_cpp = {
    # 'lcalc_Lfunction': files('lcalc_Lfunction.pyx'), # not yet on conda
}

foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/lcalc',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [],
        dependencies: [py_dep, cysignals],
    )
endforeach

