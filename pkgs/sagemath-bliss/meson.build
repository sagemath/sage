project(
    'SageMath-Bliss',
    ['c', 'cpp', 'cython'],
    version: files('../../src/VERSION.txt'),
    license: 'GPL v3',
    default_options: ['c_std=c17', 'cpp_std=c++17'],
)

# Python module
# https://mesonbuild.com/Python-module.html
py = import('python').find_installation(pure: false)

# Compilers
cc = meson.get_compiler('c')
inc_cysignals = run_command(
py,
[
    '-c',
    '''
import cysignals
print(cysignals.__file__.replace('__init__.py', ''))
    '''.strip(),
],
check: true,
).stdout().strip()
cysignals = declare_dependency(include_directories: inc_cysignals)

subdir('src')
