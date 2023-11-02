py.install_sources(
    '__init__.py',
    subdir: 'sage/libs/linbox',
)

extension_data = {
    # 'linbox_flint_interface': files('linbox_flint_interface.pyx'), # not yet on conda
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/linbox',
        install: true,
        include_directories: [],
        dependencies: [py_dep, cysignals],
    )
endforeach

