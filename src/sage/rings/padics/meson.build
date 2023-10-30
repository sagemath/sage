py.install_sources(
    'all.py',
    'eisenstein_extension_generic.py',
    'factory.py',
    'generic_nodes.py',
    'lattice_precision.py',
    'local_generic.py',
    'misc.py',
    'padic_base_generic.py',
    'padic_base_leaves.py',
    'padic_extension_generic.py',
    'padic_extension_leaves.py',
    'padic_generic.py',
    'padic_lattice_element.py',
    'padic_valuation.py',
    'precision_error.py',
    'relative_extension_leaves.py',
    'tests.py',
    'tutorial.py',
    'unramified_extension_generic.py',
    subdir: 'sage/rings/padics',
)

extension_data_cpp = {
    'padic_ZZ_pX_CA_element': files('padic_ZZ_pX_CA_element.pyx'),
    'padic_ZZ_pX_CR_element': files('padic_ZZ_pX_CR_element.pyx'),
    'padic_ZZ_pX_FM_element': files('padic_ZZ_pX_FM_element.pyx'),
    'padic_ZZ_pX_element': files('padic_ZZ_pX_element.pyx'),
    'padic_ext_element': files('padic_ext_element.pyx'),
    'pow_computer': files('pow_computer.pyx'),
    'pow_computer_ext': files('pow_computer_ext.pyx'),
}
extension_data = {
    'common_conversion': files('common_conversion.pyx'),
    'local_generic_element': files('local_generic_element.pyx'),
    'morphism': files('morphism.pyx'),
    'padic_capped_absolute_element': files('padic_capped_absolute_element.pyx'),
    'padic_capped_relative_element': files('padic_capped_relative_element.pyx'),
    'padic_fixed_mod_element': files('padic_fixed_mod_element.pyx'),
    'padic_floating_point_element': files('padic_floating_point_element.pyx'),
    'padic_generic_element': files('padic_generic_element.pyx'),
    'padic_printing': files('padic_printing.pyx'),
    #'padic_relaxed_element': files('padic_relaxed_element.pyx'), - compilation error (needs libs/linkages)
    'padic_relaxed_errors': files('padic_relaxed_errors.pyx'),
    'pow_computer_flint': files('pow_computer_flint.pyx'),
    'pow_computer_relative': files('pow_computer_relative.pyx'),
    'qadic_flint_CA': files('qadic_flint_CA.pyx'),
    'qadic_flint_CR': files('qadic_flint_CR.pyx'),
    'qadic_flint_FM': files('qadic_flint_FM.pyx'),
    'qadic_flint_FP': files('qadic_flint_FP.pyx'),
    'relative_ramified_CA': files('relative_ramified_CA.pyx'),
    'relative_ramified_CR': files('relative_ramified_CR.pyx'),
    'relative_ramified_FM': files('relative_ramified_FM.pyx'),
    'relative_ramified_FP': files('relative_ramified_FP.pyx'),
}

foreach name, pyx : extension_data
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/padics',
        install: true,
        include_directories: [inc_cpython, inc_ntl, inc_flint, inc_rings_finite, inc_ext],
        dependencies: [py_dep, cypari2, cysignals],
    )
endforeach
foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/padics',
        install: true,
        override_options : ['cython_language=cpp'],
        include_directories: [inc_cpython, inc_ntl, inc_ext],
        dependencies: [py_dep, cypari2, cysignals],
    )
endforeach
