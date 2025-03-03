py.install_sources(
  '__init__.py',
  'all.py',
  'article_heuberger_krenn_kropf_fsm-in-sage.py',
  'arxiv_0812_2725.py',
  'benchmark.py',
  'book_schilling_zabrocki_kschur_primer.py',
  'book_stein_ent.py',
  'book_stein_modform.py',
  'cmdline.py',
  'combinatorial_hopf_algebras.py',
  'finite_poset.py',
  'functools_partial_src.py',
  'gosper-sum.py',
  'lazy_imports.py',
  'modular_group_cohomology.py',
  'numpy.py',
  'parigp.py',
  'startup.py',
  'symbolic-series.py',
  'sympy.py',
  'test_deprecation.py',
  subdir: 'sage/tests',
)

extension_data = {'cython' : files('cython.pyx')}

foreach name, pyx : extension_data
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/tests',
    install: true,
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

extension_data_cpp = {'stl_vector': files('stl_vector.pyx')}

foreach name, pyx : extension_data_cpp
  py.extension_module(
    name,
    sources: pyx,
    subdir: 'sage/tests',
    install: true,
    override_options: ['cython_language=cpp'],
    include_directories: [inc_cpython, inc_rings],
    dependencies: [py_dep, cysignals, gmp],
  )
endforeach

install_subdir('books', install_dir: sage_install_dir / 'tests')
install_subdir('memcheck', install_dir: sage_install_dir / 'tests')
