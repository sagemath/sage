py.install_sources(
    'all.py',
    subdir: 'sage/rings/polynomial/weil',
)

extension_data = {
    'weil_polynomials': files('weil_polynomials.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/polynomial/weil',
        install: true,
        include_directories: [inc_flint, inc_src, inc_cpython],
        dependencies: [py_dep, cysignals],
    )
endforeach
