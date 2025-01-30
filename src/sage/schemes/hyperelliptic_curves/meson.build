inc_hypellfrob = include_directories('hypellfrob')
py.install_sources(
  '__init__.py',
  'all.py',
  'constructor.py',
  'hyperelliptic_finite_field.py',
  'hyperelliptic_g2.py',
  'hyperelliptic_generic.py',
  'hyperelliptic_padic_field.py',
  'hyperelliptic_rational_field.py',
  'invariants.py',
  'jacobian_endomorphism_utils.py',
  'jacobian_g2.py',
  'jacobian_generic.py',
  'jacobian_homset.py',
  'jacobian_morphism.py',
  'kummer_surface.py',
  'mestre.py',
  'monsky_washnitzer.py',
  subdir: 'sage/schemes/hyperelliptic_curves',
)

extension_data_cpp = {
  'hypellfrob': files(
    'hypellfrob.pyx',
    'hypellfrob/hypellfrob.cpp',
    'hypellfrob/recurrences_ntl.cpp',
  ),
}

foreach name, pyx : extension_data_cpp
  deps = [py_dep, cysignals, gmp, ntl]
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/schemes/hyperelliptic_curves',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [
      inc_cpython,
      inc_ext,
      inc_ntl,
      inc_rings,
      inc_rings_finite,
      inc_hypellfrob,
    ],
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
      install_dir: sage_install_dir / 'schemes/hyperelliptic_curves',
      install: true,
      copy: true,
    )
  endif
endforeach

