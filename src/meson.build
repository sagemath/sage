# Setup dependencies that are needed by many modules

# Meson currently ignores include_directories for Cython modules, so we
# have to add them manually.
# https://github.com/mesonbuild/meson/issues/9562
add_project_arguments('-I', meson.current_source_dir(), language : 'cython')

inc_numpy = run_command(py,
    [
    '-c',
    '''
import numpy
print(numpy.get_include())
    '''.strip()
    ],
    check: true
  ).stdout().strip()
numpy = declare_dependency(
    include_directories: inc_numpy,
)

inc_cysignals = run_command(py,
    [
    '-c',
    '''
import cysignals
print(cysignals.__file__.replace('__init__.py', ''))
    '''.strip()
    ],
    check: true
  ).stdout().strip()
cysignals = declare_dependency(
    include_directories: inc_cysignals,
)

inc_gmpy2 = run_command(py,
    [
    '-c',
    '''
import gmpy2
print(gmpy2.__file__.replace('__init__.py', ''))
    '''.strip()
    ],
    check: true
  ).stdout().strip()
gmpy2 = declare_dependency(
    include_directories: inc_gmpy2,
)

inc_cpython = include_directories('sage/cpython')
inc_rings = include_directories('sage/rings')
inc_src = include_directories('.')

# Submodules
subdir('sage/cpython')
subdir('sage/libs/gmp')
subdir('sage/misc')
subdir('sage/structure')
subdir('sage/categories')
