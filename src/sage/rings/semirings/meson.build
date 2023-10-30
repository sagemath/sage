py.install_sources(
    'all.py',
    'non_negative_integer_semiring.py',
    subdir: 'sage/rings/semirings',
)

extension_data = {
    'tropical_semiring': files('tropical_semiring.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/semirings',
        install: true,
        include_directories: [inc_cpython],
        dependencies: [py_dep],
    )
endforeach
