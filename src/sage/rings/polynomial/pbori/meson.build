py.install_sources(
    'PyPolyBoRi.py',
    '__init__.py',
    'blocks.py',
    'cnf.py',
    'easy_polynomials.py',
    'fglm.py',
    'frontend.py',
    'gbcore.py',
    'gbrefs.py',
    'heuristics.py',
    'interpolate.py',
    'interred.py',
    'll.py',
    'nf.py',
    'parallel.py',
    'randompoly.py',
    'rank.py',
    'specialsets.py',
    'statistics.py',
    subdir: 'sage/rings/polynomial/pbori',
)

extension_data_cpp = {
    'pbori': files('pbori.pyx'),
}

brial = declare_dependency(
  dependencies: [
    cc.find_library('brial', required: false),
  ]
)

foreach name, pyx : extension_data_cpp
    py.extension_module(name,
        sources: pyx,
        subdir: 'sage/rings/polynomial/pbori',
        install: true,
        include_directories: [inc_src, inc_cpython, inc_ext],
        override_options: ['cython_language=cpp'],
        dependencies: [py_dep, cysignals, brial],
    )
endforeach
