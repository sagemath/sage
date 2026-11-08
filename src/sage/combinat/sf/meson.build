py.install_sources(
  '__init__.py',
  'abreu_nigro.py',
  'all.py',
  'character.py',
  'classical.py',
  'dual.py',
  'elementary.py',
  'hall_littlewood.py',
  'hecke.py',
  'homogeneous.py',
  'jack.py',
  'k_dual.py',
  'kfpoly.py',
  'llt.py',
  'macdonald.py',
  'monomial.py',
  'multiplicative.py',
  'new_kschur.py',
  'ns_macdonald.py',
  'orthogonal.py',
  'orthotriang.py',
  'powersum.py',
  'schur.py',
  'sf.py',
  'sfa.py',
  'symplectic.py',
  'transition_kernels.pyx',
  'witt.py',
  subdir: 'sage/combinat/sf',
)

extension_data = {'transition_kernels': files('transition_kernels.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/combinat/sf',
    install: true,
    include_directories: [inc_cpython],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach
