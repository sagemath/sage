py.install_sources(
    'all.py',
    'mutable_poset.py',
    'stream.py',
    subdir: 'sage/data_structures',
)

extension_data = {
    #'binary_search': files('binary_search.pyx'), compilation error
    'bitset': files('bitset.pyx'),
    'bitset_base': files('bitset_base.pyx'),
    'blas_dict': files('blas_dict.pyx'),
    #'bounded_integer_sequences': files('bounded_integer_sequences.pyx'), requires flint
    'list_of_pairs': files('list_of_pairs.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/data_structures',
        install: true,
        include_directories: [inc_cpython],
        dependencies: [py_dep, cysignals, gmp],
    )
endforeach
