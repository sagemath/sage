py.install_sources(
  '__init__.py',
  'all.py',
  'ambient_space.py',
  'associahedron.py',
  'braid_move_calculator.py',
  'branching_rules.py',
  'cartan_matrix.py',
  'cartan_type.py',
  'coxeter_group.py',
  'coxeter_matrix.py',
  'coxeter_type.py',
  'dynkin_diagram.py',
  'extended_affine_weyl_group.py',
  'fundamental_group.py',
  'hecke_algebra_representation.py',
  'integrable_representations.py',
  'non_symmetric_macdonald_polynomials.py',
  'pieri_factors.py',
  'plot.py',
  'reflection_group_complex.py',
  'reflection_group_element.pxd',
  'reflection_group_real.py',
  'root_lattice_realization_algebras.py',
  'root_lattice_realizations.py',
  'root_space.py',
  'root_system.py',
  'type_A.py',
  'type_A_affine.py',
  'type_A_infinity.py',
  'type_B.py',
  'type_BC_affine.py',
  'type_B_affine.py',
  'type_C.py',
  'type_C_affine.py',
  'type_D.py',
  'type_D_affine.py',
  'type_E.py',
  'type_E_affine.py',
  'type_F.py',
  'type_F_affine.py',
  'type_G.py',
  'type_G_affine.py',
  'type_H.py',
  'type_I.py',
  'type_Q.py',
  'type_affine.py',
  'type_dual.py',
  'type_folded.py',
  'type_marked.py',
  'type_reducible.py',
  'type_relabel.py',
  'type_super_A.py',
  'weight_lattice_realizations.py',
  'weight_space.py',
  'weyl_characters.py',
  'weyl_group.py',
  subdir: 'sage/combinat/root_system',
)

extension_data = {
  'braid_orbit' : files('braid_orbit.pyx'),
  'reflection_group_c' : files('reflection_group_c.pyx'),
  'reflection_group_element' : files('reflection_group_element.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp]
  if name.startswith('reflection_group')
    deps += [gap]
  endif

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/root_system',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: deps,
  )

  all_deps_found = true
  foreach dep : deps
    if not dep.found()
      all_deps_found = false
      break
    endif
  endforeach
  if not all_deps_found
    configure_file(
      input: not_available,
      output: name + '.py',
      install_tag: 'python-runtime',
      install_dir: sage_install_dir / 'combinat/root_system',
      install: true,
      copy: true,
    )
  endif
endforeach

