py.install_sources(
    'all.py',
    subdir: 'sage/algebras/letterplace',
)

extension_data = {
    # 'free_algebra_element_letterplace': files('free_algebra_element_letterplace.pyx'), # not yet on conda
    # 'free_algebra_letterplace': files('free_algebra_letterplace.pyx'), # not yet on conda
    # 'letterplace_ideal': files('letterplace_ideal.pyx'), # not yet on conda
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/algebras/letterplace',
        install: true,
        include_directories: [],
        dependencies: [py_dep, cysignals],
    )
endforeach

