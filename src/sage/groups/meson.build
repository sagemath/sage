py.install_sources(
  '__init__.py',
  'all.py',
  'artin.py',
  'braid.py',
  'cactus_group.py',
  'class_function.py',
  'conjugacy_classes.py',
  'cubic_braid.py',
  'finitely_presented.py',
  'finitely_presented_catalog.py',
  'finitely_presented_named.py',
  'fqf_orthogonal.py',
  'free_group.py',
  'galois_group.py',
  'galois_group_perm.py',
  'generic.py',
  'group.pxd',
  'group_exp.py',
  'group_semidirect_product.py',
  'groups_catalog.py',
  'indexed_free_group.py',
  'kernel_subgroup.py',
  'libgap_group.py',
  'libgap_mixin.py',
  'libgap_morphism.py',
  'libgap_wrapper.pxd',
  'old.pxd',
  'pari_group.py',
  'raag.py',
  subdir: 'sage/groups',
)

extension_data = {
  'group' : files('group.pyx'),
  'libgap_wrapper' : files('libgap_wrapper.pyx'),
  'old' : files('old.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, gmp]
  if name == 'libgap_wrapper'
    deps += [gap]
  endif
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/groups',
    install: true,
    include_directories: [inc_cpython],
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
      install_dir: sage_install_dir / 'groups',
      install: true,
      copy: true,
    )
  endif
endforeach

install_subdir('abelian_gps', install_dir: sage_install_dir / 'groups')
install_subdir('additive_abelian', install_dir: sage_install_dir / 'groups')
install_subdir('affine_gps', install_dir: sage_install_dir / 'groups')
install_subdir('lie_gps', install_dir: sage_install_dir / 'groups')
subdir('matrix_gps')
install_subdir('misc_gps', install_dir: sage_install_dir / 'groups')
subdir('perm_gps')
subdir('semimonomial_transformations')
