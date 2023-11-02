py.install_sources(
    '__init__.py',
    'all.py',
    'constructor.py',
    'interface.py',
    subdir: 'sage/libs/eclib',
)

extension_data = {
    # 'homspace': files('homspace.pyx'), # not yet on conda
    # 'mat': files('mat.pyx'), # not yet on conda
    # 'mwrank': files('mwrank.pyx'), # not yet on conda
    # 'newforms': files('newforms.pyx'), # not yet on conda
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/libs/eclib',
        install: true,
        include_directories: [],
        dependencies: [py_dep, cysignals],
    )
endforeach

