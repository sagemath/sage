sage_install_dir = py.get_install_dir() / 'sage'

# Packages that need no processing and can be installed directly
no_processing = [
    'databases',
    'doctest',
    'features',
    'game_theory',
    'homology',
    'knots',
    'logic',
    'manifolds',
    'parallel',
    'repl',
    'sandpiles',
    'sat',
    'tensor',
    'topology',
    'typeset',
]
foreach package : no_processing
    install_subdir(package, install_dir: sage_install_dir)
endforeach

py.install_sources(
    '__init__.py',
    'all.py',
    'all__sagemath_bliss.py',
    'all__sagemath_categories.py',
    'all__sagemath_coxeter3.py',
    'all__sagemath_environment.py',
    'all__sagemath_mcqd.py',
    'all__sagemath_meataxe.py',
    'all__sagemath_objects.py',
    'all__sagemath_repl.py',
    'all__sagemath_sirocco.py',
    'all__sagemath_tdlib.py',
    'all_cmdline.py',
    'env.py',
    'version.py',
    subdir: 'sage',
)

subdir('cpython')
subdir('libs')
subdir('misc')
subdir('structure')
subdir('algebras')
subdir('arith')
subdir('calculus')
subdir('categories')
subdir('coding')
subdir('combinat')
subdir('crypto')
subdir('data_structures')
subdir('ext')
subdir('functions')
subdir('games')
subdir('geometry')
subdir('graphs')
subdir('groups')
subdir('interacts')
subdir('interfaces')
subdir('lfunctions')
subdir('matrix')
subdir('matroids')
subdir('modular')
subdir('modules')
subdir('monoids')
subdir('numerical')
subdir('plot')
subdir('probability')
subdir('quadratic_forms')
subdir('quivers')
subdir('rings')
subdir('schemes')
subdir('sets')
subdir('stats')
subdir('symbolic')
subdir('tests')
subdir('dynamics')
