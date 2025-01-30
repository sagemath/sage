py.install_sources(
  '__init__.py',
  'all.py',
  'category_object.pxd',
  'coerce.pxd',
  'coerce_actions.pxd',
  'coerce_dict.pxd',
  'coerce_exceptions.py',
  'coerce_maps.pxd',
  'debug_options.pxd',
  'dynamic_class.py',
  'element.pxd',
  'element_wrapper.pxd',
  'factorization.py',
  'factorization_integer.py',
  'formal_sum.py',
  'gens_py.py',
  'global_options.py',
  'indexed_generators.py',
  'list_clone.pxd',
  'list_clone_timings.py',
  'mutability.pxd',
  'nonexact.py',
  'parent.pxd',
  'parent_base.pxd',
  'parent_gens.pxd',
  'parent_old.pxd',
  'richcmp.pxd',
  'sage_object.pxd',
  'sage_object_test.py',
  'sequence.py',
  'set_factories.py',
  'set_factories_example.py',
  'support_view.py',
  'test_factory.py',
  'unique_representation.py',
  subdir: 'sage/structure',
)

extension_data = {
  'category_object' : files('category_object.pyx'),
  'coerce' : files('coerce.pyx'),
  'coerce_actions' : files('coerce_actions.pyx'),
  'coerce_dict' : files('coerce_dict.pyx'),
  'coerce_maps' : files('coerce_maps.pyx'),
  'debug_options' : files('debug_options.pyx'),
  'element' : files('element.pyx'),
  'element_wrapper' : files('element_wrapper.pyx'),
  'factory' : files('factory.pyx'),
  'list_clone' : files('list_clone.pyx'),
  'list_clone_demo' : files('list_clone_demo.pyx'),
  'list_clone_timings_cy' : files('list_clone_timings_cy.pyx'),
  'mutability' : files('mutability.pyx'),
  'parent' : files('parent.pyx'),
  'parent_base' : files('parent_base.pyx'),
  'parent_gens' : files('parent_gens.pyx'),
  'parent_old' : files('parent_old.pyx'),
  'richcmp' : files('richcmp.pyx'),
  'sage_object' : files('sage_object.pyx'),
}

foreach name, pyx : extension_data
  deps = [py_dep, cysignals, gmp]
  if name == 'coerce'
    deps += [gmpy2]
  endif

  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/structure',
    install: true,
    include_directories: [inc_cpython, inc_ext, inc_rings],
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
      install_dir: sage_install_dir / 'structure',
      install: true,
      copy: true,
    )
  endif
endforeach

install_subdir('proof', install_dir: sage_install_dir / 'structure')
