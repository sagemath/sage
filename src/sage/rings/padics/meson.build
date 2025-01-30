py.install_sources(
  '__init__.py',
  'all.py',
  'common_conversion.pxd',
  'eisenstein_extension_generic.py',
  'factory.py',
  'generic_nodes.py',
  'lattice_precision.py',
  'local_generic.py',
  'local_generic_element.pxd',
  'misc.py',
  'morphism.pxd',
  'padic_ZZ_pX_CA_element.pxd',
  'padic_ZZ_pX_CR_element.pxd',
  'padic_ZZ_pX_FM_element.pxd',
  'padic_ZZ_pX_element.pxd',
  'padic_base_generic.py',
  'padic_base_leaves.py',
  'padic_capped_absolute_element.pxd',
  'padic_capped_relative_element.pxd',
  'padic_ext_element.pxd',
  'padic_extension_generic.py',
  'padic_extension_leaves.py',
  'padic_fixed_mod_element.pxd',
  'padic_floating_point_element.pxd',
  'padic_generic.py',
  'padic_generic_element.pxd',
  'padic_lattice_element.py',
  'padic_printing.pxd',
  'padic_relaxed_element.pxd',
  'padic_relaxed_errors.pxd',
  'padic_valuation.py',
  'pow_computer.pxd',
  'pow_computer_ext.pxd',
  'pow_computer_flint.pxd',
  'pow_computer_relative.pxd',
  'precision_error.py',
  'qadic_flint_CA.pxd',
  'qadic_flint_CR.pxd',
  'qadic_flint_FM.pxd',
  'qadic_flint_FP.pxd',
  'relative_extension_leaves.py',
  'relative_ramified_CA.pxd',
  'relative_ramified_CR.pxd',
  'relative_ramified_FM.pxd',
  'relative_ramified_FP.pxd',
  'tests.py',
  'tutorial.py',
  'unramified_extension_generic.py',
  subdir: 'sage/rings/padics',
)

extension_data = {
  'common_conversion' : files('common_conversion.pyx'),
  'local_generic_element' : files('local_generic_element.pyx'),
  'morphism' : files('morphism.pyx'),
  'padic_capped_absolute_element' : files('padic_capped_absolute_element.pyx'),
  'padic_capped_relative_element' : files('padic_capped_relative_element.pyx'),
  'padic_fixed_mod_element' : files('padic_fixed_mod_element.pyx'),
  'padic_floating_point_element' : files('padic_floating_point_element.pyx'),
  'padic_generic_element' : files('padic_generic_element.pyx'),
  'padic_relaxed_element' : files('padic_relaxed_element.pyx'),
  'padic_relaxed_errors' : files('padic_relaxed_errors.pyx'),
  'qadic_flint_CA' : files('qadic_flint_CA.pyx'),
  'qadic_flint_CR' : files('qadic_flint_CR.pyx'),
  'qadic_flint_FM' : files('qadic_flint_FM.pyx'),
  'qadic_flint_FP' : files('qadic_flint_FP.pyx'),
  'relative_ramified_CA' : files('relative_ramified_CA.pyx'),
  'relative_ramified_CR' : files('relative_ramified_CR.pyx'),
  'relative_ramified_FM' : files('relative_ramified_FM.pyx'),
  'relative_ramified_FP' : files('relative_ramified_FP.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp]
  if name.startswith('qadic_flint')
    deps += [cypari2, flint]
  elif name == 'padic_capped_absolute_element' or name == 'padic_floating_point_element' or name == 'padic_fixed_mod_element' or name == 'common_conversion' or name.startswith(
    'relative_ramified',
  ) or name == 'padic_template_element' or name == 'padic_capped_relative_element'
    deps += [cypari2]
  elif name == 'padic_ext_element'
    deps += [ntl]
  elif name == 'padic_generic_element'
    deps += [
      ntl,  # Indirect dependency
    ]
  elif name == 'padic_relaxed_element'
    deps += [flint]
    if is_windows
      # error C2143: syntax error: missing ')' before '*'
      continue
    endif
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/padics',
    install: true,
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_flint,
      inc_ntl,
      inc_rings,
      inc_rings_finite,
      inc_src,
    ],
    dependencies: deps,
  )

endforeach

extension_data_cpp = {
  'padic_ZZ_pX_CA_element': files('padic_ZZ_pX_CA_element.pyx'),
  'padic_ZZ_pX_CR_element': files('padic_ZZ_pX_CR_element.pyx'),
  'padic_ZZ_pX_FM_element': files('padic_ZZ_pX_FM_element.pyx'),
  'padic_ZZ_pX_element': files('padic_ZZ_pX_element.pyx'),
  'padic_ext_element': files('padic_ext_element.pyx'),
  'padic_printing': files('padic_printing.pyx'),
  'pow_computer': files('pow_computer.pyx'),
  'pow_computer_ext': files('pow_computer_ext.pyx'),
  'pow_computer_flint': files('pow_computer_flint.pyx'),
  'pow_computer_relative': files('pow_computer_relative.pyx'),
}

foreach name, pyx : extension_data_cpp
  deps = [py_dep, cysignals, gmp, ntl]
  if name == 'pow_computer_flint'
    deps += [flint]
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/rings/padics',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_flint,
      inc_ntl,
      inc_rings,
      inc_rings_finite,
    ],
    dependencies: deps,
  )
  )
endforeach

