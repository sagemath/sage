py.install_sources(
  'all.py',
  'conway_polynomials.py',
  'element_base.pxd',
  'element_givaro.pxd',
  'element_ntl_gf2e.pxd',
  'element_pari_ffelt.pxd',
  'finite_field_base.pxd',
  'finite_field_constructor.py',
  'finite_field_givaro.py',
  'finite_field_ntl_gf2e.py',
  'finite_field_pari_ffelt.py',
  'finite_field_prime_modn.py',
  'galois_group.py',
  'hom_finite_field.pxd',
  'hom_finite_field_givaro.pxd',
  'hom_prime_finite_field.pxd',
  'homset.py',
  'integer_mod.pxd',
  'integer_mod_limits.h',
  'integer_mod_ring.py',
  'maps_finite_field.py',
  'residue_field.pxd',
  'stdint.pxd',
  subdir: 'sage/rings/finite_rings',
)

extension_data = {
  'element_base' : files('element_base.pyx'),
  'element_pari_ffelt' : files('element_pari_ffelt.pyx'),
  'finite_field_base' : files('finite_field_base.pyx'),
  'hom_finite_field' : files('hom_finite_field.pyx'),
  'hom_prime_finite_field' : files('hom_prime_finite_field.pyx'),
  'integer_mod' : files('integer_mod.pyx'),
  'residue_field' : files('residue_field.pyx'),
  'residue_field_givaro' : files('residue_field_givaro.pyx'),
  'residue_field_ntl_gf2e' : files('residue_field_ntl_gf2e.pyx'),
  'residue_field_pari_ffelt' : files('residue_field_pari_ffelt.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp]
  if name == 'element_pari_ffelt'
    deps += [
      cypari2,
      pari,  # Runtime dependency
    ]
  elif name == 'residue_field_pari_ffelt'
    deps += [cypari2]
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/finite_rings',
    install: true,
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_ntl,
      inc_rings,
      inc_rings_finite,
    ],
    dependencies: deps,
  )
endforeach

extension_data_cpp = {
  'element_givaro': files('element_givaro.pyx'),
  'element_ntl_gf2e': files('element_ntl_gf2e.pyx'),
  'hom_finite_field_givaro': files('hom_finite_field_givaro.pyx'),
}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/finite_rings',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_ntl,
      inc_rings,
      inc_rings_finite,
    ],
    dependencies: [py_dep, cypari2, cysignals, givaro, gmp, m, ntl, pari],
  )
endforeach

