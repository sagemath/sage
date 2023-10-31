sage_install_dir = py.get_install_dir() / 'sage'

# Packages that need no processing and can be installed directly
no_processing = [
    'repl',
    'doctest',
    'features',
    'tensor',
    'topology',
    'typeset',
]
foreach package : no_processing
    install_subdir(package, install_dir: sage_install_dir)
endforeach

py.install_sources(
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
subdir('categories')
subdir('sets')
subdir('data_structures')
subdir('rings')
subdir('arith')
subdir('combinat')
subdir('modules')
subdir('tests')
subdir('ext')
subdir('interfaces')
