py.install_sources(
    'all.py',
    'conway_polynomials.py',
    'finite_field_constructor.py',
    'finite_field_givaro.py',
    'finite_field_ntl_gf2e.py',
    'finite_field_pari_ffelt.py',
    'finite_field_prime_modn.py',
    'galois_group.py',
    'homset.py',
    'integer_mod_ring.py',
    'maps_finite_field.py',
    subdir: 'sage/rings/finite_rings',
)

extension_data_cpp = {
    'element_givaro': files('element_givaro.pyx'),
    'element_ntl_gf2e': files('element_ntl_gf2e.pyx'),
    'hom_finite_field': files('hom_finite_field.pyx'),
    'hom_finite_field_givaro': files('hom_finite_field_givaro.pyx'),
}
extension_data = {
    'element_base': files('element_base.pyx'),
    'element_pari_ffelt': files('element_pari_ffelt.pyx'),
    'finite_field_base': files('finite_field_base.pyx'),
    'hom_prime_finite_field': files('hom_prime_finite_field.pyx'),
    'integer_mod': files('integer_mod.pyx'),
    'residue_field': files('residue_field.pyx'),
    'residue_field_givaro': files('residue_field_givaro.pyx'),
    'residue_field_ntl_gf2e': files('residue_field_ntl_gf2e.pyx'),
    'residue_field_pari_ffelt': files('residue_field_pari_ffelt.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/finite_rings',
        install: true,
        include_directories: [inc_cpython, inc_ntl, inc_rings],
        dependencies: [py_dep, cysignals, cypari2],
    )
endforeach
foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/finite_rings',
        install: true,
        override_options: ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_ntl, inc_rings, inc_ext],
        dependencies: [py_dep, cysignals, cypari2],
    )
endforeach
