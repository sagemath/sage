py.install_sources(
  'all.py',
  'buzzard.py',
  'congroup.py',
  'congroup_element.py',
  'cusps.py',
  'cusps_nf.py',
  'dims.py',
  'dirichlet.py',
  'etaproducts.py',
  'hypergeometric_motive.py',
  'multiple_zeta.py',
  'multiple_zeta_F_algebra.py',
  subdir: 'sage/modular',
)

extension_data = {'hypergeometric_misc' : files('hypergeometric_misc.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/modular',
    install: true,
    include_directories: [inc_cpython],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

install_subdir('abvar', install_dir: sage_install_dir / 'modular')
subdir('arithgroup')
install_subdir('btquotients', install_dir: sage_install_dir / 'modular')
install_subdir('drinfeld_modform', install_dir: sage_install_dir / 'modular')
install_subdir('hecke', install_dir: sage_install_dir / 'modular')
install_subdir('local_comp', install_dir: sage_install_dir / 'modular')
subdir('modform')
install_subdir(
  'modform_hecketriangle',
  install_dir: sage_install_dir / 'modular',
)
subdir('modsym')
install_subdir('overconvergent', install_dir: sage_install_dir / 'modular')
subdir('pollack_stevens')
install_subdir('quasimodform', install_dir: sage_install_dir / 'modular')
install_subdir('quatalg', install_dir: sage_install_dir / 'modular')
install_subdir('ssmod', install_dir: sage_install_dir / 'modular')
