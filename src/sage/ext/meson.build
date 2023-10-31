py.install_sources(
    'all__sagemath_objects.py',
    subdir: 'sage/ext',
)

extension_data = {
    'fast_callable': files('fast_callable.pyx'),
    'fast_eval': files('fast_eval.pyx'),
    'memory': files('memory.pyx'),
    'memory_allocator': files('memory_allocator.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/ext',
        install: true,
        include_directories: [inc_cpython],
        dependencies: [py_dep, cysignals, cypari2, gmp],
    )
endforeach

subdir('interpreters')
