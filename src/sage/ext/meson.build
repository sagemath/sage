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

interpreters = custom_target(
    'sage.ext.interpreters',
    output: 'all.py',
    input: '../../sage_setup/autogen/interpreters/__init__.py',
    command: [py, '-m', 'sage_setup.autogen.interpreters', meson.current_source_dir() / 'interpreters'],
    # This is actually against the Meson philosophy, which does not like in-source builds.
    # So normally you would do something like the following (left for a follow-up)
    # command: [py, '-m', 'sage_setup.autogen.interpreters','@OUTDIR@'],
    env: ['PYTHONPATH=' + meson.current_source_dir() / '..' / '..'],
    install: true,
    install_dir: 'interpreters'
)

subdir('interpreters')
