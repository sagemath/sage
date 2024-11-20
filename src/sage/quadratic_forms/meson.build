py.install_sources(
  'all.py',
  'binary_qf.py',
  'bqf_class_group.py',
  'constructions.py',
  'extras.py',
  'qfsolve.py',
  'quadratic_form.py',
  'quadratic_form__automorphisms.py',
  'quadratic_form__count_local_2.py',
  'quadratic_form__equivalence_testing.py',
  'quadratic_form__genus.py',
  'quadratic_form__local_density_congruence.py',
  'quadratic_form__local_density_interfaces.py',
  'quadratic_form__local_field_invariants.py',
  'quadratic_form__local_normal_form.py',
  'quadratic_form__local_representation_conditions.py',
  'quadratic_form__mass.py',
  'quadratic_form__mass__Conway_Sloane_masses.py',
  'quadratic_form__mass__Siegel_densities.py',
  'quadratic_form__neighbors.py',
  'quadratic_form__reduction_theory.py',
  'quadratic_form__siegel_product.py',
  'quadratic_form__split_local_covering.py',
  'quadratic_form__ternary_Tornaria.py',
  'quadratic_form__theta.py',
  'quadratic_form__variable_substitutions.py',
  'random_quadraticform.py',
  'special_values.py',
  'ternary_qf.py',
  subdir: 'sage/quadratic_forms',
)

extension_data = {
  'count_local_2' : files('count_local_2.pyx'),
  'quadratic_form__evaluate' : files('quadratic_form__evaluate.pyx'),
  'ternary' : files('ternary.pyx'),
}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/quadratic_forms',
    install: true,
    include_directories: [],
    dependencies: [py_dep, gmp],
  )
endforeach

install_subdir('genera', install_dir: sage_install_dir / 'quadratic_forms')
